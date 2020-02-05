import React from "react"
import { Col, FormGroup, Input, Label } from 'reactstrap';
import { Field } from 'formik';
import { FieldErrorMessage, ModelCardNew } from '../../../../genui';
import * as Yup from 'yup';

function DrugExNetValidationFields(props) {
  const validationStrategyPrefix = props.validationStrategyPrefix;

  return (
    <React.Fragment>
      <FormGroup row>
        <Label htmlFor={`${validationStrategyPrefix}.validSetSize`} sm={4}>Validation Set Size</Label>
        <Col sm={8}>
          <Field name={`${validationStrategyPrefix}.validSetSize`} as={Input} type="number" step="1"/>
        </Col>
      </FormGroup>
      <FieldErrorMessage name={`${validationStrategyPrefix}.validSetSize`}/>
    </React.Fragment>
  )
}

function DrugExNetExtraFields(props) {
  return (
    <React.Fragment>
      <FormGroup>
        <Label htmlFor="molset">Corpus</Label>
        <Field name="molset" as={Input} type="select">
          {
            props.molsets.map((molset) => <option key={molset.id} value={molset.id}>{molset.name}</option>)
          }
        </Field>
      </FormGroup>
      <FieldErrorMessage name="molset"/>

      <FormGroup>
        <Label htmlFor="parent">Parent Network</Label>
        <Field name="parent" as={Input} type="select">
          {
            props.models.map((model) => <option key={model.id} value={model.id}>{model.name}</option>)
          }
        </Field>
      </FormGroup>
      <FieldErrorMessage name="parent"/>
    </React.Fragment>
  )
}

export class DrugExNetCreateCard extends React.Component {

  render() {
    let molsets = [];
    Object.keys(this.props.compoundSets).forEach(
      (key) => molsets = molsets.concat(this.props.compoundSets[key])
    );

    const validationStrategyInit = {
      validSetSize: 0,
    };
    const extraParamInit = {
      parent: this.props.models.length > 0 ? this.props.models[0] : undefined,
      molset: molsets[0].id,
    };

    const validationStrategySchema = {
      validSetSize: Yup.number().integer().min(0, 'Validation set size must be positive or zero.'),
    };
    const extraParamsSchema = {
      parent: Yup.number().integer().positive('Parent network ID must be a positive integer.'),
      molset: Yup.number().integer().positive('Molecule set ID must be a positive integer.').required('You need to supply a set of compounds to create the corpus from.')
    };

    return (
      <ModelCardNew
        {...this.props}
        molsets={molsets}
        validationStrategyInit={validationStrategyInit}
        extraParamsInit={extraParamInit}
        validationStrategySchema={validationStrategySchema}
        extraParamsSchema={extraParamsSchema}
        validationStrategyFields={DrugExNetValidationFields}
        extraFields={DrugExNetExtraFields}
      />
    )
  }
}

export class DrugExAgentCreateCard extends React.Component {

  render() {
    return <div>Agents Card New</div>;
  }
}