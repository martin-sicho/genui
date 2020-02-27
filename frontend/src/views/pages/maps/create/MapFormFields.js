import React from 'react';
import { Col, FormGroup, Input, Label } from 'reactstrap';
import { Field } from 'formik';
import { FieldErrorMessage } from '../../../../genui';

export function MapTrainFields(props) {
  const trainingStrategyPrefix = props.trainingStrategyPrefix;

  return (
    <React.Fragment>
      <FormGroup>
        <Label htmlFor={`${trainingStrategyPrefix}.descriptors`}>Descriptor Sets</Label>
        <p>
          Choose one or more descriptor sets to use in the calculations.
        </p>
        <Field name={`${trainingStrategyPrefix}.descriptors`} as={Input} type="select" multiple>
          {
            props.descriptors.map((desc) => (
              <option key={desc.id} value={desc.id}>
                {desc.name}
              </option>
            ))
          }
        </Field>
      </FormGroup>
      <FieldErrorMessage name={`${trainingStrategyPrefix}.descriptors`}/>
    </React.Fragment>
  )
}

export function MapExtraFields(props) {
  const molsets = props.molsets;

  return (
    <React.Fragment>
      <FormGroup row>
        <Label htmlFor="molsets" sm={4}>Compound Sets</Label>
        <Col sm={8}>
          <Field name="molsets" as={Input} type="select" multiple>
            {
              molsets.map(molset => (
                <option key={molset.id} value={molset.id}>
                  {molset.name}
                </option>
              ))
            }
          </Field>
        </Col>
      </FormGroup>
      <FieldErrorMessage name="molsets"/>
    </React.Fragment>
  )
}