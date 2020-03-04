import React from 'react';
import * as Yup from 'yup';
import { CardBody, CardHeader, FormGroup, Input, Label } from 'reactstrap';
import { Field } from 'formik';
import { FieldErrorMessage, ComponentWithObjects, GenericNewMolSetCard } from '../../../../genui';

function ExtraFormFields(props) {
  return (
    <React.Fragment>

      <FormGroup>
        <Label htmlFor="source">Generator</Label>
        <Field name="source" as={Input} type="select">
          {
            props.generators.map((source) => <option key={source.id} value={source.id}>{source.name}</option>)
          }
        </Field>
        <FieldErrorMessage name="source"/>
      </FormGroup>

      <FormGroup>
        <Label htmlFor="nSamples">Number of compounds to generate</Label>
        <Field name="nSamples" as={Input} type="number"/>
      </FormGroup>
      <FieldErrorMessage name="nSamples"/>
    </React.Fragment>
  )
}

function NoGeneratorsCard(props) {
  return (
    <React.Fragment>
      <CardHeader>No Generators Available</CardHeader>
      <CardBody>
        There are currently no generators available. You have to create one first.
      </CardBody>
    </React.Fragment>
  )
}

function GeneratedCardNew(props) {

  const extraInitVals = {
    source : undefined,
    nSamples : 100
  };

  const extraValidSchemas = {
    source: Yup.number().integer().positive("Source generator ID has to be a positive number.").required('Source generator ID is required.'),
    nSamples: Yup.number().min(1, 'Required number of generated compounds must be at least 1.').required('You must provide the number of samples.')
  };

  return (
    <ComponentWithObjects
      {...props}
      objectListURL={new URL('all/', props.apiUrls.generatorsRoot)}
      emptyClassName="Generator"
      render={
        (generators) => {
          let flattened = [];
          Object.keys(generators).forEach(key => {
            flattened = flattened.concat(generators[key]);
          });
          extraInitVals.source = flattened.length > 0 ? flattened[0].id : undefined;
          return (
            flattened.length > 0 ? <GenericNewMolSetCard
              {...props}
              generators={flattened}
              cardHeader="Generate Compounds"
              extraFormInitVals={extraInitVals}
              extraFormValidSchemas={extraValidSchemas}
              additionalFieldsComponent={ExtraFormFields}
            /> : <NoGeneratorsCard/>
          )
        }
      }
    />
  )
}

export default GeneratedCardNew;