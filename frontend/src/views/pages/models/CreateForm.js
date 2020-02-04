import React from "react"
import { Button, CardBody, CardFooter, Col, Form, FormGroup, Input, Label } from 'reactstrap';
import { Field, Formik } from 'formik';
import * as Yup from 'yup';
import { FieldErrorMessage } from '../../../genui';

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

class ModelForm extends React.Component {

  constructor(props) {
    super(props);

    this.state = {
      formIsSubmitting : false,
    }
  }

  setFormSubmitting = (state) => {
    this.setState({formIsSubmitting : state});
  };

  render() {
    const formIsSubmitting = this.state.formIsSubmitting;
    const validationSchema = this.props.validationSchema;
    const initialValues = this.props.initialValues;
    const validationStrategyPrefix = "validationStrategy";
    const trainingStrategyPrefix = "trainingStrategy";

    return (
      <React.Fragment>
        <CardBody className="scrollable">
          <Formik
            initialValues={initialValues}
            validationSchema={validationSchema}
            onSubmit={
              (values) => {
                this.setFormSubmitting(true);
                this.props.handleCreate(values);
              }
            }
          >
            {
              formik => (
                <Form id="model-create-form" onSubmit={formik.handleSubmit} className="unDraggable">
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
                    <Label htmlFor="molset">Training Set</Label>
                    <Field name="molset" as={Input} type="select">
                      {
                        this.props.molsets.map((molset) => <option key={molset.id} value={molset.id}>{molset.name}</option>)
                      }
                    </Field>
                  </FormGroup>
                  <FieldErrorMessage name="molset"/>

                  <h4>Training Parameters</h4>

                  <FormGroup>
                    <Label htmlFor={`${trainingStrategyPrefix}.mode`}>Mode</Label>
                    <Field name={`${trainingStrategyPrefix}.mode`} as={Input} type="select">
                      {
                        this.props.modes.map((mode) => <option key={mode.id} value={mode.id}>{mode.name}</option>)
                      }
                    </Field>
                    <FieldErrorMessage name={`${trainingStrategyPrefix}.mode`}/>
                  </FormGroup>

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
                        this.props.descriptors.map((desc) => <option key={desc.id} value={desc.id}>{desc.name}</option>)
                      }
                    </Field>
                  </FormGroup>
                  <FieldErrorMessage name={`${trainingStrategyPrefix}.descriptors`}/>

                  {/*{this.props.parameters.length > 0 ? <h4>Algorithm Parameters</h4> : null}*/}

                  {
                    this.props.parameters.map(param => {
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

                  {formik.initialValues.hasOwnProperty(validationStrategyPrefix) > 0 ?
                    <React.Fragment>
                      <h4>Validation Parameters</h4>

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

                      <FormGroup row>
                        <Label htmlFor={`${validationStrategyPrefix}.metrics`} sm={4}>Validation Metrics</Label>
                        <Col sm={8}>
                          <Field name={`${validationStrategyPrefix}.metrics`} as={Input} type="select" multiple>
                            {
                              this.props.metrics.map(metric => (
                                <option key={metric.id} value={metric.id}>
                                  {metric.name}
                                </option>
                              ))
                            }
                          </Field>
                        </Col>
                      </FormGroup>
                      <FieldErrorMessage name={`${validationStrategyPrefix}.metrics`}/>
                    </React.Fragment>
                    : null
                  }

                  <FormGroup>

                  </FormGroup>
                </Form>
              )
            }
          </Formik>
        </CardBody>
        <CardFooter>
          <Button block form="model-create-form" type="submit" color="primary" disabled={formIsSubmitting}>{formIsSubmitting ? "Creating..." : "Create"}</Button>
        </CardFooter>
      </React.Fragment>
    )
  }
}

class ModelCreateForm extends React.Component {
  CTYPE_TO_VALIDATOR = {
    string: Yup.string().required(),
    integer: Yup.number().required(),
    float: Yup.number().required(),
    bool: Yup.bool().required()
  };

  CTYPE_TO_DEFAULT = {
    string: "",
    integer: 1,
    float: 1.0,
    bool: false
  };

  constructor(props) {
    super(props);

    this.chosenAlgorithm = this.props.chosenAlgorithm;
    this.modes = this.props.chosenAlgorithm.validModes;
    this.parameters = this.props.chosenAlgorithm.parameters;
    this.molsets = this.props.molsets;
    this.metrics = this.props.metrics;
    this.initialValues = this.generateInit();
    this.schema = this.generateSchema();
  }

  generateInit = () => {
    const trainingStrategyDefaultInit = {
      algorithm: this.chosenAlgorithm.id,
      mode: this.modes[0].id,
    };
    const trainingStrategyInit = Object.assign(trainingStrategyDefaultInit, this.props.trainingStrategyInit);
    const validationStrategyDefaultInit = {
      metrics: [this.metrics[0].id]
    };
    const validationStrategyInit = Object.assign(validationStrategyDefaultInit, this.props.validationStrategyInit);

    let initialValues = {
      name: `New ${this.chosenAlgorithm.name} Model`,
      description: '',
      molset: this.molsets[0].id,
      project: this.props.project.id,
      trainingStrategy: trainingStrategyInit,
      validationStrategy: validationStrategyInit
    };
    initialValues = Object.assign(initialValues, this.props.extraParamsInit);
    const parameterDefaults = {};
    for (const param of this.parameters) {
      parameterDefaults[param.name] = this.CTYPE_TO_DEFAULT[param.contentType]
    }
    initialValues.trainingStrategy.parameters = parameterDefaults;

    // console.log(initialValues);

    return initialValues;
  };

  generateSchema = () => {
    const parameterValidators = {};
    for (const param of this.parameters) {
      parameterValidators[param.name] = this.CTYPE_TO_VALIDATOR[param.contentType]
    }
    const trainingStrategyDefault = {
      algorithm: Yup.number().integer().positive("Algorithm ID needs to be a positive number").required('Algorithm ID must be supplied'),
      mode: Yup.number().integer()
        .max(256, 'Mode must be 256 characters or less.').required('You must specify a mode.'),
      parameters: Yup.object().shape(parameterValidators)
    };
    const trainingStrategy = Object.assign(trainingStrategyDefault, this.props.trainingStrategySchema);

    const validationStrategyDefault = {
      metrics: Yup.array().of(Yup.number().positive('Metric ID must be a positive integer.')).required('You need to supply one or more validation metrics for validation.'),
    };
    const validationStrategy = Object.assign(validationStrategyDefault, this.props.validationStrategySchema);

    let validationObj = {
      name: Yup.string()
        .max(256, 'Name must be less than 256 characters long.')
        .required('Name is required.'),
      description: Yup.string()
        .max(10000, 'Description must be 10,000 characters or less.'),
      project: Yup.number().integer().positive("Project ID needs to be a positive number").required('Project ID must be supplied'),
      trainingStrategy: Yup.object().shape(
        trainingStrategy
      ),
      validationStrategy: Yup.object().shape(
        validationStrategy
      )
    };
    validationObj = Object.assign(validationObj, this.props.extraParamsSchema);
    // console.log(validationObj);
    return Yup.object().shape(validationObj);
  };

  render() {
    return (
      <ModelForm
        {...this.props}
        initialValues={this.initialValues}
        validationSchema={this.schema}
        modes={this.modes}
        parameters={this.parameters}
      />
    );
  }
}

export default ModelCreateForm;