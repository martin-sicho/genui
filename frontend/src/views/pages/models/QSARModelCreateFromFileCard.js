import React from "react";
import { FormikModelUploadForm, ModelCardNew } from '../../../genui';
import * as Yup from 'yup';
import { DescriptorsField, PredictionsFields } from './QSARModelFormFields';

export default function QSARModelCreateFromFileCard (props) {
  const trainingStrategyInit = {
    descriptors: [props.descriptors[0].id],
  };

  const extraParamsInit = {
    predictionsType: "",
    predictionsUnits: ""
  };

  const trainingStrategySchema = {
    descriptors: Yup.array().of(Yup.number().positive('Descriptor set ID must be a positive integer.')).required('You need to supply one or more descriptor sets for training.'),
  };

  const extraParamsSchema = {
    predictionsType: Yup.string().required('Predictions activity type has to be set.').min(1, "Predictions type cannot be empty.").max(128, 'Predictions type name cannot be longer than 128 characters.'),
    predictionsUnits: Yup.string().max(128, 'Predictions units name cannot be longer than 128 characters.')
  };

  return (
    <ModelCardNew
      {...props}
      form={FormikModelUploadForm}
      formNameSuffix="create-upload"
      omitAlgParams={true}
      omitValidation={true}
      enableFileUploads={true}
      extraParamsInit={extraParamsInit}
      extraParamsSchema={extraParamsSchema}
      trainingStrategyInit={trainingStrategyInit}
      trainingStrategySchema={trainingStrategySchema}
      trainingStrategyFields={(props) => (
        <DescriptorsField
          {...props}
          description="Check the descriptor sets that were used to train this model."
        />
      )}
      extraFields={PredictionsFields}
      onValuesInit={(values, state) => {
        if (state.modes) {
          const mode = state.modes[0];
          if (mode.name === "classification") {
            values.predictionsType = "Active Probability"
          }
        }
        return values;
      }}
      prePost={(data) => {
        if (data.predictionsUnits === "") {
          data.predictionsUnits = null;
        }
        return data;
      }}
    />)
}