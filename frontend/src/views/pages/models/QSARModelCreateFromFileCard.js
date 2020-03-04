import React from "react";
import { FormikModelUploadForm, ModelCardNew } from '../../../genui';

export default function QSARModelCreateFromFileCard (props) {
  return (
    <ModelCardNew
      {...props}
      form={FormikModelUploadForm}
      formNameSuffix="create-upload"
      omitAlgParams={true}
      omitValidation={true}
      enableFileUploads={true}
    />)
}