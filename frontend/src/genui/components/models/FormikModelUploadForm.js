import React from 'react';
import { Field, Formik } from 'formik';
import { Form, FormGroup, FormText, Input, Label } from 'reactstrap';
import { FieldErrorMessage, FileUpload } from '../../index';

export  default function FormikModelUploadForm (props) {
  const trainingStrategyPrefix = "trainingStrategy";
  const modes = props.modes;

  const TrainingStrategyExtras = props.trainingStrategyFields;
  const ExtraFields = props.extraFields;
  return (
    <Formik
      initialValues={props.initialValues}
      validationSchema={props.validationSchema}
      onSubmit={props.onSubmit}
    >
      {
        formik => (
          <Form id={`${props.modelClass}-${props.formNameSuffix}-form`} onSubmit={formik.handleSubmit} className="unDraggable">
            <FormGroup>
              <Label htmlFor="name">Model Name</Label>
              <Field name="name" as={Input} type="text"/>
            </FormGroup>
            <FieldErrorMessage name="name"/>
            <FormGroup>
              <Label htmlFor="description">Description</Label>
              <Field name="description" as={Input} type="textarea" placeholder="Write more about this model if needed..."/>
            </FormGroup>
            <FieldErrorMessage name="description"/>
            <Field name={`${trainingStrategyPrefix}.algorithm`} as={Input} type="number" hidden/>
            <Field name="project" as={Input} type="number" hidden/>

            <FormGroup>
              <Label htmlFor="modelFile">Model File</Label>
              <Field
                name="modelFile"
                component={FileUpload}
              />
              <FormText color="muted">
                Upload a model file. If you upload a model file, the model will not be trained, but rather an
                attempt will be made to load it using the information below. Make sure to use a supported
                file extension in the uploaded file name. It will be used to infer a proper deserialization
                procedure.
              </FormText>
            </FormGroup>
            <FieldErrorMessage
              name="modelFile"
            />

            <h4>Model Info</h4>

            {ExtraFields ? <ExtraFields {...props}/> : null}

            <FormGroup>
              <Label htmlFor={`${trainingStrategyPrefix}.mode`}>Mode</Label>
              <Field name={`${trainingStrategyPrefix}.mode`} as={Input} type="select">
                {
                  modes.map((mode) => <option key={mode.id} value={mode.id}>{mode.name}</option>)
                }
              </Field>
              <FieldErrorMessage name={`${trainingStrategyPrefix}.mode`}/>
            </FormGroup>

            {
              TrainingStrategyExtras ?
                <TrainingStrategyExtras
                  {...props}
                  trainingStrategyPrefix={trainingStrategyPrefix}
                /> : null
            }
          </Form>
        )
      }
    </Formik>
  )
}