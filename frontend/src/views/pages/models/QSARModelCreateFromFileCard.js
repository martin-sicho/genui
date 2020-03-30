import React from "react";
import { FormikModelUploadForm, ModelCardNew } from '../../../genui';
import * as Yup from 'yup';
import { DescriptorsField } from './QSARModelFormFields';

export default function QSARModelCreateFromFileCard (props) {
  const trainingStrategyInit = {
    descriptors: [props.descriptors[0].id],
  };

  const trainingStrategySchema = {
    descriptors: Yup.array().of(Yup.number().positive('Descriptor set ID must be a positive integer.')).required('You need to supply one or more descriptor sets for training.'),
  };

  return (
    <ModelCardNew
      {...props}
      form={FormikModelUploadForm}
      formNameSuffix="create-upload"
      omitAlgParams={true}
      omitValidation={true}
      enableFileUploads={true}
      trainingStrategyInit={trainingStrategyInit}
      trainingStrategySchema={trainingStrategySchema}
      trainingStrategyFields={(props) => (
        <DescriptorsField
          {...props}
          description="Check the descriptor sets that were used to train this model."
        />
      )}
    />)
}