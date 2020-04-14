import React from 'react';
import { Col, FormGroup, Input, Label } from 'reactstrap';
import { Field } from 'formik';
import { FieldErrorMessage } from '../../../genui';

export function DescriptorsField (props) {
  const trainingStrategyPrefix = props.trainingStrategyPrefix;
  const description = props.description;

  return (
    <React.Fragment>
      <FormGroup>
        <Label htmlFor={`${trainingStrategyPrefix}.descriptors`}>Descriptor Sets</Label>
        <p>
          {description}
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

export function PredictionsFields(props) {
  return (
    <React.Fragment>
      <FormGroup>
        <Label htmlFor="predictionsType">Activity Type for Predictions</Label>
        <p>
          This will be the activity type for the output values of this QSAR model.
          It is recommended to leave this at the default value. If you are uploading a model,
          this value should be set to the activity type value of the uploaded model.
        </p>
        <Field name="predictionsType" as={Input} type="text"/>
      </FormGroup>
      <FieldErrorMessage name="predictionsType"/>

      <FormGroup>
        <Label htmlFor="predictionsUnits">Activity Units for Predictions</Label>
        <p>
          Use this to specify the dimension of the model output activity type.
          Leave it blank to determine this automatically or if the activity type has no dimension.
        </p>
        <Field name="predictionsUnits" as={Input} type="text"/>
      </FormGroup>
      <FieldErrorMessage name="predictionsUnits"/>
    </React.Fragment>
  )
}

export function QSARTrainingFields (props) {
  const trainingStrategyPrefix = props.trainingStrategyPrefix;

  return (
    <React.Fragment>
      <FormGroup>
        <Label htmlFor={`${trainingStrategyPrefix}.activitySet`}>Activity Set</Label>
        <Field name={`${trainingStrategyPrefix}.activitySet`} as={Input} type="select">
          {
            props.activitySets.map((set) => <option key={set.id} value={set.id}>{set.name}</option>)
          }
        </Field>
      </FormGroup>
      <FieldErrorMessage name={`${trainingStrategyPrefix}.activitySet`}/>

      <FormGroup>
        <Label htmlFor={`${trainingStrategyPrefix}.activityType`}>Activity Type</Label>
        <Field name={`${trainingStrategyPrefix}.activityType`} as={Input} type="select">
          {
            props.activityTypes.map(type => <option key={type.id} value={type.id}>{type.value}</option>)
          }
        </Field>
      </FormGroup>
      <FieldErrorMessage name={`${trainingStrategyPrefix}.activityType`}/>

      {
        props.modes.find(mode => mode.name === "classification") ? (
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
          </React.Fragment>
        ) : null
      }

      <DescriptorsField {...props} description="Choose one or more descriptor sets to use in the calculations."/>
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
        <Label htmlFor="molset">Compound Set</Label>
        <p>Compounds from this set and their bioactivity data will be used for training.</p>
        <Field name="molset" as={Input} type="select">
          {
            molsets.map((molset) => <option key={molset.id} value={molset.id}>{molset.name}</option>)
          }
        </Field>
      </FormGroup>
      <FieldErrorMessage name="molset"/>

      <PredictionsFields {...props}/>
    </React.Fragment>
  )
}