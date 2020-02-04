import React from "react";
import {
  ComponentWithObjects,
  ComponentWithResources,
  ModelCardNew,
  ModelsPage,
  FieldErrorMessage,
} from '../../../genui';
import ModelCard from './ModelCard';
import * as Yup from 'yup';
import { Col, FormGroup, Input, Label } from 'reactstrap';
import { Field } from 'formik';

function QSARTrainingFields (props) {
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
        <Field name={`${trainingStrategyPrefix}.activityThreshold`} as={Input} type="number"/>
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

function QSARValidationFields(props) {
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

function QSARExtraFields(props) {
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

function QSARModelCreateCard (props) {
  let molsets = [];
  Object.keys(props.compoundSets).forEach(
    (key) => molsets = molsets.concat(props.compoundSets[key])
  );

  const trainingStrategyInit = {
    activityThreshold : 6.5,
    descriptors: [props.descriptors[0].id],
  };
  const validationStrategyInit = {
    cvFolds: 10,
    validSetSize: 0.2,
  };
  const extraParamInit = {
    molset: molsets[0].id,
  };

  const trainingStrategySchema = {
    activityThreshold: Yup.number().min(0, 'Activity threshold must be zero or positive.'),
    descriptors: Yup.array().of(Yup.number().positive('Descriptor set ID must be a positive integer.')).required('You need to supply one or more descriptor sets for training.'),
  };
  const validationStrategySchema = {
    cvFolds: Yup.number().integer().min(0, 'Number of CV folds must be at least 0.'),
    validSetSize: Yup.number().min(0.0, 'Validation set size must be at least 0.0.').max(1.0,'Validation set size is expressed as a fraction, which needs to be less than 1.0.'),
  };

  const extraParamsSchema = {
    molset: Yup.number().integer().positive('Molecule set ID must be a positive integer.').required('You need to supply a training set of compounds.'),
  };

  return (
    <ModelCardNew
      {...props}
      molsets={molsets}
      trainingStrategyInit={trainingStrategyInit}
      validationStrategyInit={validationStrategyInit}
      extraParamsInit={extraParamInit}
      trainingStrategySchema={trainingStrategySchema}
      validationStrategySchema={validationStrategySchema}
      extraParamsSchema={extraParamsSchema}
      trainingStrategyFields={QSARTrainingFields}
      validationStrategyFields={QSARValidationFields}
      extraFields={QSARExtraFields}
    />)
}

function Models(props) {
  const resources = {
    algorithmChoices : new URL('algorithms/', props.apiUrls.qsarRoot),
    descriptors: new URL('descriptors/', props.apiUrls.qsarRoot),
    metrics: new URL('metrics/', props.apiUrls.qsarRoot)
  };
  return (
    <ComponentWithResources definition={resources}>
      {
        (allLoaded, resources) => (
          allLoaded ? <ComponentWithObjects
            objectListURL={new URL('all/', props.apiUrls.compoundSetsRoot)}
            {...props}
            render={
              (
                ...args
              ) => {
                const [compoundSets] = [...args];
                const compoundSetsAvailable = !(Object.keys(compoundSets).length === 0 && compoundSets.constructor === Object);
                return (compoundSetsAvailable ? <ModelsPage
                  {...props}
                  {...resources}
                  modelClass="QSARModel"
                  listURL={new URL(`models/`, props.apiUrls.qsarRoot)}
                  modelComponent={ModelCard} // TODO: rename this and put some common parts in the genui package
                  newModelComponent={QSARModelCreateCard}
                  compoundSets={compoundSets}
                /> : <div><p>There are currently no compound sets. You need to create one before building a QSAR model.</p></div>)
              }
            }
          /> : <div>Loading...</div>
        )
      }
    </ComponentWithResources>
  );
}

export default Models;