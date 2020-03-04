import React from "react";
import { Field, Formik } from 'formik';
import { Button, CardBody, CardFooter, Form, FormGroup, Input, Label } from 'reactstrap';
import { FieldErrorMessage } from '../../index';
import * as Yup from 'yup';

class GenericNewMolSetForm extends React.Component {
  render() {
    const AdditionalFields = this.props.additionalFieldsComponent;
    return (
      <Formik
        initialValues={this.props.initialValues}
        validationSchema={this.props.validationSchema}
        onSubmit={this.props.onSubmit}
      >
        {
          formik => (
            <Form id={this.props.formID} onSubmit={formik.handleSubmit} className="unDraggable">
              <FormGroup>
                <Label htmlFor="name">Name</Label>
                <Field name="name" as={Input} type="text"/>
              </FormGroup>
              <FieldErrorMessage name="name"/>
              <FormGroup>
                <Label htmlFor="description">Description</Label>
                <Field name="description" as={Input} type="textarea"/>
              </FormGroup>
              <FieldErrorMessage name="description"/>

              <AdditionalFields {...this.props} formik={formik}/>
            </Form>
          )
        }
      </Formik>
    )
  }
}

export function NewMolSetFormRenderer(props) {
  const [formIsSubmitting, setFormIsSubmitting] = React.useState(false);

  const initialValues = Object.assign({
    name: 'New Compound Set Name',
    description: 'Detailed description of the compound set...'
  }, props.extraFormInitVals);
  const validationSchema = Yup.object(
    Object.assign({
      name: Yup.string()
        .max(256, 'Name must be less than 256 character long.')
        .required('Name is required.'),
      description: Yup.string()
        .max(10000, 'Description must be 10,000 characters or less'),
    }, props.extraFormValidSchemas)
  );

  const id = `${props.currentMolsetClass}-create-form`;
  return (
    <React.Fragment>
      <CardBody className="scrollable">
        <GenericNewMolSetForm
          {...props}
          formID={id}
          initialValues={initialValues}
          validationSchema={validationSchema}
          onSubmit={
            (values) => {
              setFormIsSubmitting(true);
              props.handleCreate(values);
            }
          }
        />
      </CardBody>
      <CardFooter>
        <Button block form={id} type="submit" color="primary" disabled={formIsSubmitting}>{formIsSubmitting ? "Creating..." : "Create"}</Button>
      </CardFooter>
    </React.Fragment>
  );
}

export default NewMolSetFormRenderer;