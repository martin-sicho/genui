import * as Yup from 'yup';
import { ModelCardNew } from '../../../genui';
import React from 'react';
import { QSARExtraFields, QSARTrainingFields, QSARValidationFields } from './QSARModelFormFields';

export default function QSARModelCreateCard (props) {
  let molsets = [];
  Object.keys(props.compoundSets).forEach(
    (key) => molsets = molsets.concat(props.compoundSets[key])
  );

  const trainingStrategyInit = {
    activityThreshold : 6.5,
    descriptors: [props.descriptors[0].id],
  };
  const validationStrategyInit = {
    cvFolds: 10,
    validSetSize: 0.2,
  };
  const extraParamInit = {
    molset: molsets[0].id,
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
    <ModelCardNew
      {...props}
      molsets={molsets}
      trainingStrategyInit={trainingStrategyInit}
      validationStrategyInit={validationStrategyInit}
      extraParamsInit={extraParamInit}
      trainingStrategySchema={trainingStrategySchema}
      validationStrategySchema={validationStrategySchema}
      extraParamsSchema={extraParamsSchema}
      trainingStrategyFields={QSARTrainingFields}
      validationStrategyFields={QSARValidationFields}
      extraFields={QSARExtraFields}
    />)
}