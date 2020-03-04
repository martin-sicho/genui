import React from "react"
import { Field } from 'formik';
import { UncontrolledAlert } from 'reactstrap';

const FieldErrorMessage = ({ name }) => (
  <Field
    name={name}
  >
    {
      ({field, form, meta}) => {
        const error = meta.error;
        const touch = meta.touched;
        if (touch && error) {
          return <UncontrolledAlert color="danger">{error}</UncontrolledAlert>;
        }

        return null;
      }
    }
  </Field>
);

export default FieldErrorMessage;