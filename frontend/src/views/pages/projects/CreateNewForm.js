import { useFormik } from 'formik';
import * as Yup from 'yup';
import { Button, Form, FormGroup, Input, Label, UncontrolledAlert } from 'reactstrap';
import React from 'react';

export function CreateNewForm(props) {
  const formik = useFormik({
    initialValues: {
      name: 'New Project Name',
      description: 'Short description of the project...',
    },
    validationSchema: Yup.object({
      name: Yup.string()
        .max(256, 'Name must be less than 256 character long.')
        .required('Name is required.'),
      description: Yup.string()
        .max(10000, 'Description must be 10,000 characters or less'),
    }),
    onSubmit: values => {
      props.handleCreate(values);
    },
  });

  return (
    <Form onSubmit={formik.handleSubmit} className="unDraggable">
      <FormGroup>
        <Label htmlFor="name">Name</Label>
        <Input
          {...formik.getFieldProps('name')}
          type="text"
        />
      </FormGroup>
      {formik.touched.name && formik.errors.name ?
        <UncontrolledAlert color="danger">{formik.errors.name}</UncontrolledAlert> : null}
      <FormGroup>
        <Label htmlFor="description">Description</Label>
        <Input
          {...formik.getFieldProps('description')}
          type="textarea"
        />
      </FormGroup>
      {formik.touched.description && formik.errors.description ?
        <UncontrolledAlert color="danger">{formik.errors.description}</UncontrolledAlert> : null}
      <Button type="submit" color="primary">Create & Open</Button>
    </Form>
  );
}