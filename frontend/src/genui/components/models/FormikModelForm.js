import React from 'react';
import { Field, Formik } from 'formik';
import { Col, Form, FormGroup, Input, Label } from 'reactstrap';
import { FieldErrorMessage } from '../../index';

class ParameterField extends React.Component {
  CTYPE_TO_COMPONENT = {
    string: name => <Field name={name} as={Input} type="text"/>,
    integer: name => <Field name={name} as={Input} type="number"/>,
    float: name => <Field name={name} as={Input} type="number" step="0.01"/>,
    bool: name => <Field name={name} as={Input} type="checkbox"/>
  };

  render() {
    return this.CTYPE_TO_COMPONENT[this.props.parameter.contentType](this.props.name)
  }
}

function FormikModelForm (props) {
  const validationStrategyPrefix = "validationStrategy";
  const trainingStrategyPrefix = "trainingStrategy";
  const modes = props.modes;
  const parameters = props.parameters;
  const metrics = props.metrics;

  const TrainingStrategyExtras = props.trainingStrategyFields;
  const ValidationStrategyExtras = props.validationStrategyFields;
  const ExtraFields = props.extraFields;
  return (
    <Formik
      initialValues={props.initialValues}
      validationSchema={props.validationSchema}
      onSubmit={props.onSubmit}
    >
      {
        formik => (
          <Form id={`${props.modelClass}-create-form`} onSubmit={formik.handleSubmit} className="unDraggable">
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

            {ExtraFields ? <ExtraFields {...props}/> : null}

            <h4>Training Parameters</h4>

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

            {parameters.length > 0 ? <h4>{props.chosenAlgorithm.name} Parameters</h4> : null}

            {
              parameters.map(param => {
                const name = `${trainingStrategyPrefix}.parameters.${param.name}`;
                return (
                  <FormGroup key={name} row>
                    <Label htmlFor={name} sm={4}>{param.name}</Label>
                    <Col sm={8}>
                      <ParameterField parameter={param} name={name}/>
                      <FieldErrorMessage name={name}/>
                    </Col>
                  </FormGroup>
                )})
            }

            {formik.initialValues.hasOwnProperty(validationStrategyPrefix) ?
              <React.Fragment>
                <h4>Validation Parameters</h4>

                <FormGroup row>
                  <Label htmlFor={`${validationStrategyPrefix}.metrics`} sm={4}>Validation Metrics</Label>
                  <Col sm={8}>
                    <Field name={`${validationStrategyPrefix}.metrics`} as={Input} type="select" multiple>
                      {
                        metrics.map(metric => (
                          <option key={metric.id} value={metric.id}>
                            {metric.name}
                          </option>
                        ))
                      }
                    </Field>
                  </Col>
                </FormGroup>
                <FieldErrorMessage name={`${validationStrategyPrefix}.metrics`}/>

                {
                  ValidationStrategyExtras ?
                  <ValidationStrategyExtras
                    validationStrategyPrefix={validationStrategyPrefix}
                  /> : null
                }
              </React.Fragment>
              : null
            }

            <FormGroup>

            </FormGroup>
          </Form>
        )
      }
    </Formik>
  )
}

export default FormikModelForm;