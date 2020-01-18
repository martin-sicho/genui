import React from "react"
import { Button, CardBody, CardFooter, Form, FormGroup, Input, Label} from 'reactstrap';
import { Field, Formik } from 'formik';
import * as Yup from 'yup';
import { FieldErrorMessage } from '../../../genui';

export function ModelCreateForm(props) {

  const [formIsSubmitting, setFormIsSubmitting] = React.useState(false);

  return (
    <React.Fragment>
      <CardBody className="scrollable">
        <Formik
          initialValues={{
            name: 'New Model',
            description: 'Write more about this model if needed...',
          }}
          validationSchema={Yup.object(
            {
              name: Yup.string()
                .max(256, 'Name must be less than 256 characters long.')
                .required('Name is required.'),
              description: Yup.string()
                .max(10000, 'Description must be 10,000 characters or less'),
            })
          }
          onSubmit={
            (values) => {
              setFormIsSubmitting(true);
              props.handleCreate(values);
            }
          }
        >
          {
            formik => (
              <Form id="model-create-form" onSubmit={formik.handleSubmit} className="unDraggable">
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

              {/*  TODO: add fields*/}

              </Form>
            )
          }
        </Formik>
      </CardBody>
      <CardFooter>
        <Button block form="model-create-form" type="submit" color="primary" disabled={formIsSubmitting}>{formIsSubmitting ? "Creating..." : "Create"}</Button>
      </CardFooter>
    </React.Fragment>
  );
}