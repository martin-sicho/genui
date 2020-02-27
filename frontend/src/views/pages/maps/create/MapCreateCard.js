import React from 'react';
import { ModelCardNew } from '../../../../genui';
import * as Yup from 'yup';
import { MapExtraFields, MapTrainFields } from './MapFormFields';

export default function MapCreateCard(props) {
  let molsets = [];
  Object.keys(props.compoundSets).forEach(
    (key) => molsets = molsets.concat(props.compoundSets[key])
  );

  const trainingStrategyInit = {
    descriptors: [props.descriptors[0].id],
  };
  const extraParamInit = {
    molsets: [molsets[0].id],
  };

  const trainingStrategySchema = {
    activityThreshold: Yup.number().min(0, 'Activity threshold must be zero or positive.'),
    descriptors: Yup.array().of(Yup.number().positive('Descriptor set ID must be a positive integer.')).required('You need to supply one or more descriptor sets for training.'),
  };

  const extraParamsSchema = {
    molsets: Yup.array().of(
      Yup.number().positive(
        'Molecule set ID must be a positive integer.'
      )).required('You have to select at least one compound set to map.'),
  };

  return (
    <ModelCardNew
      {...props}
      molsets={molsets}
      trainingStrategyInit={trainingStrategyInit}
      extraParamsInit={extraParamInit}
      trainingStrategySchema={trainingStrategySchema}
      extraParamsSchema={extraParamsSchema}
      trainingStrategyFields={MapTrainFields}
      extraFields={MapExtraFields}
    />)
}