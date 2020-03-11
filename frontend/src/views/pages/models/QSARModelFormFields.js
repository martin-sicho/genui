import React from 'react';
import { Col, FormGroup, Input, Label } from 'reactstrap';
import { Field } from 'formik';
import { FieldErrorMessage } from '../../../genui';

export function QSARTrainingFields (props) {
  const trainingStrategyPrefix = props.trainingStrategyPrefix;

  return (
    <React.Fragment>
      <FormGroup>
        <Label htmlFor={`${trainingStrategyPrefix}.activityThreshold`}>Activity Threshold</Label>
        <p>
          This is only relevant in classification mode.
          Molecules with their primary activity measure
          higher than or equal to this value will be considered active.
        </p>
        <Field name={`${trainingStrategyPrefix}.activityThreshold`} as={Input} type="number" step={0.01}/>
      </FormGroup>
      <FieldErrorMessage name={`${trainingStrategyPrefix}.activityThreshold`}/>

      <FormGroup>
        <Label htmlFor={`${trainingStrategyPrefix}.descriptors`}>Descriptor Sets</Label>
        <p>
          Choose one or more descriptor sets to use in the calculations.
        </p>
        <Field name={`${trainingStrategyPrefix}.descriptors`} as={Input} type="select" multiple>
          {
            props.descriptors.map((desc) => <option key={desc.id} value={desc.id}>{desc.name}</option>)
          }
        </Field>
      </FormGroup>
      <FieldErrorMessage name={`${trainingStrategyPrefix}.descriptors`}/>
    </React.Fragment>
  )
}

export function QSARValidationFields(props) {
  const validationStrategyPrefix = props.validationStrategyPrefix;

  return (
    <React.Fragment>
      <FormGroup row>
        <Label htmlFor={`${validationStrategyPrefix}.cvFolds`} sm={4}>Cross-Validation Folds</Label>
        <Col sm={8}>
          <Field name={`${validationStrategyPrefix}.cvFolds`} as={Input} type="number"/>
        </Col>
      </FormGroup>
      <FieldErrorMessage name={`${validationStrategyPrefix}.cvFolds`}/>

      <FormGroup row>
        <Label htmlFor={`${validationStrategyPrefix}.validSetSize`} sm={4}>Validation Set Size</Label>
        <Col sm={8}>
          <Field name={`${validationStrategyPrefix}.validSetSize`} as={Input} type="number" step="0.01"/>
        </Col>
      </FormGroup>
      <FieldErrorMessage name={`${validationStrategyPrefix}.validSetSize`}/>
    </React.Fragment>
  )
}

export function QSARExtraFields(props) {
  const molsets = props.molsets;

  return (
    <React.Fragment>
      <FormGroup>
        <Label htmlFor="molset">Training Set</Label>
        <Field name="molset" as={Input} type="select">
          {
            molsets.map((molset) => <option key={molset.id} value={molset.id}>{molset.name}</option>)
          }
        </Field>
      </FormGroup>
      <FieldErrorMessage name="molset"/>
    </React.Fragment>
  )
}