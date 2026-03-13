import json
import importlib
import tempfile
from abc import ABC

from rdkit import Chem
import numpy as np
import pandas as pd

from qsprpred.data import QSPRDataset
from qsprpred import models as qsprpred_models
from qsprpred.models import CrossValAssessor, TestSetAssessor
from sklearn.model_selection import KFold, StratifiedKFold

from django.core.exceptions import ImproperlyConfigured
from genui.compounds.models import Molecule

from .qsprpred_utils import RecordProgressMonitor, MetricsAggregator, build_split_instance
from genui.compounds.models import ActivityTypes, ActivitySet
from genui.models import models as core_models
from genui.models.genuimodels.bases import Algorithm, PredictionMixIn, ValidationMixIn, ProgressMixIn, ModelBuilder
from genui.qsar import models as qsar_models
from .bases import EmbeddingBuilderMixIn

clustering_module = importlib.import_module("qsprpred.data.chem.clustering")


class BasicQSARModelBuilder(EmbeddingBuilderMixIn, PredictionMixIn, ValidationMixIn, ProgressMixIn, ModelBuilder, ABC):
    def __init__(self, instance: qsar_models.Model, progress=None, onFitCall=None, validations=None):
        super().__init__(instance, progress, onFitCall)
        self.validations = validations if validations and len(
            validations) > 0 else self.instance.trainingStrategy.validationStrategies.all()
        self.temp_dir = tempfile.TemporaryDirectory()
        self.dataset = None

    def build(self) -> qsar_models.QSARModel:
        if not self.validations:
            raise ImproperlyConfigured("You cannot build a QSAR model without validation strategies.")

        if not self.molsets:
            raise ImproperlyConfigured("You cannot build a QSAR model without an associated molecule set.")

        self.progressStages = [
            "Fetching activities...",
            "Calculating embeddings..."
        ]
        if self.hyper_param_opt:
            self.progressStages.append("Optimizing hyperparameters...")
            if self.hyper_param_opt.__class__.__name__ == "GridSearchOptimization":
                num_combinations = 1
                for params in self.search_space.values():
                    num_combinations *= len(params)
                for i in range(num_combinations):
                    self.progressStages.append(f"Evaluating combination {i}...")
            elif self.hyper_param_opt.__class__.__name__ == "OptunaOptimization":
                for i in range(self.hyper_param_opt.nTrials):
                    self.progressStages.append(f"Evaluating trial {i}...")

        for validation in self.validations:
            self.progressStages.extend([f"CV fold {x + 1}" for x in range(validation.cvFolds)])
        self.progressStages.extend(["Fitting model on the training set...", "Validating on test set..."])
        self.progressStages.extend(["Fitting the final model...", "Saving the model..."])

        self.recordProgress()
        dataset = self.saveActivities()

        self.recordProgress()
        dataset.prepareDataset(feature_calculators=self.embeddingCalculators)

        monitor = RecordProgressMonitor(self.recordProgress)
        monitor.on_fold_start = True
        monitor.on_optimization_start = True
        metrics_aggregator = MetricsAggregator(self.metricFunctions[0][0], [m for m in self.metricFunctions], self, monitor,
                                               core_models.ModelPerformanceCV)
        self._model = self.algorithmClass(self)

        for validation_strategy_index, validation in enumerate(self.validations):
            monitor.validation_index = validation_strategy_index
            if hasattr(validation, 'dataSplit') and hasattr(validation, 'cvFolds'):
                split_instance = build_split_instance(validation.dataSplit, dataset=dataset)
                dataset.split(split_instance)

                if self.hyper_param_opt:
                    monitor.on_assessment_start = True
                    monitor.on_fold_start = False
                    metrics_aggregator.perfClass = core_models.HyperparameterOptimizationPerformance
                    optimizer_class = getattr(qsprpred_models, self.hyper_param_opt.__class__.__name__)
                    kwargs = {"param_grid": self.search_space,
                              "model_assessor": CrossValAssessor(self.hypo_metric),
                              "score_aggregation": self.hyper_param_aggregator,
                              "monitor": monitor}
                    if "nTrials" in dir(self.hyper_param_opt):
                        kwargs["n_trials"] = self.hyper_param_opt.nTrials
                    optimizer = optimizer_class(**kwargs)
                    params = optimizer.optimize(self.model.model, dataset)
                    old_params = json.loads(self.model.params['parameters'])
                    old_params.update(**params)
                    self.model.params["parameters"] = json.dumps(old_params)
                    metrics_aggregator.perfClass = core_models.ModelPerformanceCV
                    monitor.on_fold_start = True
                    monitor.on_assessment_start = False

                # Perform cross-validation
                is_regression = self.training.mode.name == Algorithm.REGRESSION
                if is_regression:
                    folds = KFold(validation.cvFolds)
                else:
                    folds = StratifiedKFold(validation.cvFolds)
                cva = CrossValAssessor(metrics_aggregator,
                                       folds,
                                       monitor,
                                       use_proba=self.training.mode.name == Algorithm.CLASSIFICATION, )
                cva(self.model.model, dataset, False, )

        monitor.on_assessment_start = True
        monitor.validation_index = -1 # Test set performance
        metrics_aggregator.perfClass = core_models.ModelPerformance
        tsa = TestSetAssessor(metrics_aggregator,
                              monitor,
                              use_proba=self.training.mode.name == Algorithm.CLASSIFICATION, )
        tsa(self.model.model, dataset, False)

        self.recordProgress()
        self.model.model.fitDataset(dataset, save_model=False)

        # Final validation (optional, as it's not truly a validation)
        # y_predicted = self.model.predict(dataset.X)
        # for validation in self.validations:
        # self.validate(dataset.y, y_predicted)  # Not a validation

        self.recordProgress()
        self.saveFile()
        return self.instance

    def populateActivitySet(self, aset: qsar_models.ModelActivitySet):
        if not self.instance.predictionsType:
            raise Exception("The activity type for QSAR model predictions is not specified.")

        aset.activities.all().delete()
        molecules = aset.molecules.molecules.all()
        predictions = self.predictMolecules(molecules)

        for mol, prediction in zip(molecules, predictions):
            qsar_models.ModelActivity.objects.create(
                value=prediction,
                type=self.instance.predictionsType,
                units=self.instance.predictionsUnits,
                source=aset,
                molecule=mol,
            )

        return aset.activities.all()

    def predictMolecules(self, mols):
        smiles = []
        failed_indices = []
        for idx, mol in enumerate(mols):
            if mol:
                smiles.append(Chem.MolToSmiles(mol) if type(mol) == Chem.Mol else mol)
            else:
                failed_indices.append(idx)

        predictions = [-1] * len(mols)
        if len(failed_indices) == len(mols):
            return np.array(predictions)

        real_predictions = list(self.predictMols(smiles))
        for idx, prediction in enumerate(predictions):
            if idx not in failed_indices:
                predictions[idx] = real_predictions.pop(0)
        assert len(real_predictions) == 0
        return np.array(predictions)

    def getDataset(self) -> QSPRDataset:
        return self.dataset

    def saveActivities(self):
        if not self.getDataset():
            activity_set = ActivitySet.objects.get(pk=self.training.activitySet.id)
            activity_type = self.training.activityType
            if not activity_set:
                raise Exception("No activity set specified.")
            if not activity_type:
                raise Exception("No activity type specified.")

            compounds, activities, units = activity_set.cleanForModelling(activity_type)
            if not len(compounds) == len(activities):
                raise Exception(
                    f'Number of compounds in a QSAR model ({len(compounds)}) is different from the set of activities assigned to them ({len(activities)}). Something went wrong when the data was cleaned for modeling.')

            if self.training.mode.name == Algorithm.CLASSIFICATION:
                if not self.instance.predictionsType:
                    self.instance.predictionsType = ActivityTypes.objects.get_or_create(
                        value="Active Probability"
                    )[0]

            if not self.instance.predictionsType:
                self.instance.predictionsType = activity_type
            if not self.instance.predictionsUnits:
                self.instance.predictionsUnits = units

            self.instance.save()
            if self.training.mode.name == Algorithm.CLASSIFICATION:
                target_props = {"name": activity_type.value, "task": "SINGLECLASS",
                                "th": [self.training.activityThreshold], }
            elif self.training.mode.name == Algorithm.MULTICLASS:
                target_props = {"name": activity_type.value, "task": "MULTICLASS", "th": self.training.activityThreshold}
            else:
                target_props = {"name": activity_type.value, "task": "REGRESSION"}
            smiles = self.getSmiles(compounds)
            df = pd.DataFrame({"SMILES": smiles, activity_type.value: activities})
            self.dataset = QSPRDataset(self.instance.name, [target_props], df, store_dir=self.temp_dir.name)
            return self.dataset

    def load_model(self):  # Backup, for future endeavours
        self._model = self.model.load_model(self.instance.files.filter(note=self.model.model_name)[0].path)

    def getSmiles(self, mols=None) -> list:
        if mols is not None:
            smiles = [x.canonicalSMILES if isinstance(x, Molecule) else x for x in mols]
        elif hasattr(self, "mols"):
            smiles = [x.canonicalSMILES if isinstance(x, Molecule) else x for x in self.mols]
        elif hasattr(self, "molsets"):
            smiles = []
            for molset in self.molsets:
                for mol in molset.molecules.all():
                    smiles.append(mol.canonicalSMILES)
        else:
            raise Exception("No molecules to calculate embeddings from.")

        return smiles

    def predictMols(self, smiles=None):
        smiles = self.getSmiles(smiles)
        return self.model.predictMols(smiles)
