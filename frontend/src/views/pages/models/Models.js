import React from "react";
import { ComponentWithObjects, ComponentWithResources, ModelCardNew, ModelsPage, ModelFormRenderer } from '../../../genui';
import ModelCard from './ModelCard';
import ModelForm from './CreateForm';
import * as Yup from 'yup';

function QSARModelForm (props) {
  const trainingStrategyInit = {
    activityThreshold : 6.5,
    descriptors: [props.descriptors[0].id],
  };
  const validationStrategyInit = {
    cvFolds: 10,
    validSetSize: 0.2,
  };
  const extraParamInit = {
    molset: props.molsets[0].id,
  };

  const trainingStrategySchema = {
    activityThreshold: Yup.number().min(0, 'Activity threshold must be zero or positive.'),
    descriptors: Yup.array().of(Yup.number().positive('Descriptor set ID must be a positive integer.')).required('You need to supply one or more descriptor sets for training.'),
  };
  const validationStrategySchema = {
    cvFolds: Yup.number().integer().min(0, 'Number of CV folds must be at least 0.'),
    validSetSize: Yup.number().min(0.0, 'Validation set size must be at least 0.0.').max(1.0,'Validation set size is expressed as a fraction, which needs to be less than 1.0.'),
  };

  const extraParamsSchema = {
    molset: Yup.number().integer().positive('Molecule set ID must be a positive integer.').required('You need to supply a training set of compounds.'),
  };

  return (
    <ModelFormRenderer
      {...props}
      trainingStrategyInit={trainingStrategyInit}
      validationStrategyInit={validationStrategyInit}
      extraParamsInit={extraParamInit}
      trainingStrategySchema={trainingStrategySchema}
      validationStrategySchema={validationStrategySchema}
      extraParamsSchema={extraParamsSchema}
      formComponent={ModelForm}
    />)
}

function NewQSARCard (props) {
  let molsets = [];
  Object.keys(props.compoundSets).forEach(
    (key) => molsets = molsets.concat(props.compoundSets[key])
  );

  return <ModelCardNew {...props} molsets={molsets} formComponent={QSARModelForm}/>
}

function Models(props) {
  const resources = {
    algorithmChoices : new URL('algorithms/', props.apiUrls.qsarRoot),
    descriptors: new URL('descriptors/', props.apiUrls.qsarRoot),
    metrics: new URL('metrics/', props.apiUrls.qsarRoot)
  };
  return (
    <ComponentWithResources definition={resources}>
      {
        (allLoaded, resources) => (
          allLoaded ? <ComponentWithObjects
            objectListURL={new URL('all/', props.apiUrls.compoundSetsRoot)}
            {...props}
            render={
              (
                ...args
              ) => {
                const [compoundSets] = [...args];
                const compoundSetsAvailable = !(Object.keys(compoundSets).length === 0 && compoundSets.constructor === Object);
                return (compoundSetsAvailable ? <ModelsPage
                  {...props}
                  {...resources}
                  modelClass="QSARModel"
                  listURL={new URL(`models/`, props.apiUrls.qsarRoot)}
                  modelComponent={ModelCard}
                  newModelComponent={NewQSARCard}
                  compoundSets={compoundSets}
                /> : <div><p>There are currently no compound sets. You need to create one before building a QSAR model.</p></div>)
              }
            }
          /> : <div>Loading...</div>
        )
      }
    </ComponentWithResources>
  );
}

export default Models;