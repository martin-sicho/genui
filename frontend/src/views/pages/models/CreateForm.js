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
                    <Field name="description" as={Input} type="textarea"/>
                  </FormGroup>
                  <FieldErrorMessage name="description"/>

                  {formik.initialValues.hasOwnProperty("mode") ?
                    <FormGroup>
                      <Label htmlFor="mode">Mode</Label>
                      <Field name="mode" as={Input} type="select">
                        {
                          this.props.modes.map((mode) => <option key={mode.name} value={mode.name}>{mode.name}</option>)
                        }
                      </Field>
                      <FieldErrorMessage name="mode"/>
                    </FormGroup>
                    : null
                  }

                  <FormGroup>
                    <Label htmlFor="activityThrs">Activity Threshold</Label>
                    <p>
                      This is only relevant in classification mode.
                      Molecules with their primary activity measure
                      higher than or equal to this value will be considered active.
                    </p>
                    <Field name="activityThrs" as={Input} type="number"/>
                  </FormGroup>
                  <FieldErrorMessage name="activityThrs"/>

                  {this.props.parameters.length > 0 ? <h4>Algorithm Parameters</h4> : null}

                  {
                    this.props.parameters.map(param => {
                      const name = `parameters.${param.name}`;
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

                  {/*TODO: add descriptors and validation parameters*/}
                  {/*{this.props.validationParams.length > 0 ? <h4>Training Parameters</h4> : null}*/}

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
    string: "Some String",
    integer: 0,
    float: 0.0,
    bool: false
  };

  constructor(props) {
    super(props);

    this.newModel = this.props.newModel;
    this.modes = this.props.newModel.validModes;
    this.parameters = this.props.newModel.parameters;
  }

  render() {
    let initialValues = {
      name: `New ${this.newModel.name} Model`,
      description: 'Write more about this model if needed...',
      mode: this.modes[0].name,
      activityThrs : 6.5,
      // molset: molsets[0].id
    };
    const parameterDefaults = {parameters : {}};
    for (const param of this.parameters) {
      parameterDefaults.parameters[param.name] = this.CTYPE_TO_DEFAULT[param.contentType]
    }
    initialValues = Object.assign(
      initialValues,
      parameterDefaults
    );

    let validationObj = {
      name: Yup.string()
        .max(256, 'Name must be less than 256 characters long.')
        .required('Name is required.'),
      description: Yup.string()
        .max(10000, 'Description must be 10,000 characters or less.'),
      mode: Yup.string()
        .max(256, 'Mode must be 256 characters or less.'),
      activityThrs: Yup.number().min(0, 'Activity threshold must be zero or positive.'),
      // molset: Yup.number().min(1, 'Molecule set ID must be a positive integer.')
    };
    const parameterValidators = {};
    for (const param of this.parameters) {
      parameterValidators[param.name] = this.CTYPE_TO_VALIDATOR[param.contentType]
    }
    validationObj = Object.assign(
      validationObj,
      {parameters : Yup.object().shape(parameterValidators)}
    );
    const validationSchema = Yup.object().shape(validationObj);

    return (
      <ModelForm
        initialValues={initialValues}
        validationSchema={validationSchema}
        modes={this.modes}
        parameters={this.parameters}
        molsets={this.props.molsets}
        modelDefinition={this.newModel}
        handleCreate={this.props.handleCreate}
      />
    );
  }
}

export default ModelCreateForm;