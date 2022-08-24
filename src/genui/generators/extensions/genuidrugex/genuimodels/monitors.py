"""
monitors

Created by: Martin Sicho
On: 30-01-20, 15:05
"""
import torch
from django.core.files.base import ContentFile
from drugex.training.monitors import DictMonitor

from genui.generators.extensions.genuidrugex.genuimodels.metrics import DrugExLoss, SMILESErrorRate, SMILESUniqueRate, \
    MeanDrExDesirability
from genui.generators.extensions.genuidrugex.models import ModelPerformanceDrugEx
from genui.models.models import ModelFile


class DrugExMonitor(DictMonitor):

    def __init__(self, model_instance, progress_callback):
        super().__init__(verbose=True, clean_after_epoch=True)
        self.modelInstance = model_instance
        self.progressCall = progress_callback
        self.moleculesFile = ModelFile.create(self.modelInstance, 'molecules.tsv', ContentFile(''), kind=ModelFile.AUXILIARY, note='molecules')
        self.headerWritten = False

    def saveMolecules(self, df):
        open_mode = 'a' if self.headerWritten else 'w'
        df.to_csv(
            self.moleculesFile.path,
            sep='\t',
            index=False,
            header=not self.headerWritten,
            mode=open_mode,
            encoding='utf-8',
            na_rep='NA'
        )
        self.headerWritten = True

    def close(self):
        pass

    def saveModel(self, model):
        super(DrugExMonitor, self).saveModel(model)
        model_file = self.modelInstance.modelFile
        if not model_file:
            model_file = ModelFile.create(self.modelInstance, 'best.pkg', ContentFile(''), kind=ModelFile.MAIN, note='main')
        torch.save(self.bestState, model_file.path)

    def saveEpochData(self, df):
        loss_train = df['mean_train_loss'][0]
        loss_valid = df['loss_valid'][0]
        error_rate = 1 - df['valid_ratio'][0]
        unique_ratio = df['unique_ratio'][0]
        desire_ratio = df['desire_ratio'][0]
        if loss_train is not None:
            self.savePerformance(DrugExLoss.getDjangoModel(), loss_train, False)
        if loss_valid is not None:
            self.savePerformance(DrugExLoss.getDjangoModel(), loss_valid, True)
        if error_rate is not None:
            self.savePerformance(SMILESErrorRate.getDjangoModel(), error_rate, True)
        if unique_ratio:
            self.savePerformance(SMILESUniqueRate.getDjangoModel(), unique_ratio, True)
        if desire_ratio:
            self.savePerformance(MeanDrExDesirability.getDjangoModel(), desire_ratio, True)
        self.progressCall(self.currentEpoch-1)

    def savePerformance(self, metric, value, isValidation, note=""):
        return ModelPerformanceDrugEx.objects.create(
            metric=metric,
            value=value,
            isOnValidationSet=isValidation,
            model=self.modelInstance,
            epoch=self.currentEpoch-1,
            step=self.currentEpoch-1,
            note=note
        )