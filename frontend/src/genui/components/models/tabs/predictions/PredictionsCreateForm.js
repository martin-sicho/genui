import React from 'react';
import { Button, Form, FormGroup, Input, Label } from 'reactstrap';
import { Field, Formik } from 'formik';
import { FieldErrorMessage } from '../../../../index';
import * as Yup from 'yup';

function PredictionsCreateForm(props) {
  const molsets = props.molsets;
  const [isSubmitting, setIsSubmitting] = React.useState(false);
  return (
    <Formik
      initialValues={{
        name: '',
        description: '',
        molecules: molsets.length > 0 ? molsets[0].id : undefined
      }}
      validationSchema={Yup.object({
        name: Yup.string()
          .max(256, 'Name must be less than 256 character long.')
          .required('Name is required.'),
        description: Yup.string()
          .max(10000, 'Description must be 10,000 characters or less'),
        molecules: Yup.number().integer().positive(
          'Molecule set ID must be a positive integer.'
        ).required(
          'You need to supply a set of compounds to make predictions for.'
        ),
      })}
      onSubmit={
        (values) => {
          setIsSubmitting(true);
          fetch(
            props.predictionsListUrl
            , {
              method: 'POST'
              , body: JSON.stringify(values)
              , headers: {
                'Content-Type': 'application/json'
              },
              credentials: "include",
            }
          ).then((data) => props.handleResponseErrors(data, "Creating predictions failed. Data wrong or incomplete?"))
            .then(data => {
              props.handleAddActivitySet(props.defaultClass, data);
              setIsSubmitting(false);
            })
            .catch(
            error => console.log(error)
          );
        }
      }
    >
      {
        formik => (
          <Form onSubmit={formik.handleSubmit} className="unDraggable">
            <FormGroup>
              <Label htmlFor="name">Name</Label>
              <Field name="name" as={Input} type="text" placeholder="Name for the set of activities (will be shown in summaries and tables)"/>
            </FormGroup>
            <FieldErrorMessage name="name"/>
            <FormGroup>
              <Label htmlFor="description">Description</Label>
              <Field name="description" as={Input} type="textarea" placeholder="Notes about the new set of activities..."/>
            </FormGroup>
            <FieldErrorMessage name="description"/>

            <FormGroup>
              <Label htmlFor="molecules">Compound Set to Predict</Label>
              <Field name="molecules" as={Input} type="select">
                {
                  molsets.map((molset) => <option key={molset.id} value={molset.id}>{molset.name}</option>)
                }
              </Field>
            </FormGroup>
            <FieldErrorMessage name="molecules"/>

            <Button block type="submit" color="primary" disabled={isSubmitting}>{isSubmitting ? "Creating..." : "Create"}</Button>
          </Form>
        )
      }
    </Formik>
  )
}

export default PredictionsCreateForm;