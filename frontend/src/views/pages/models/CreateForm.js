import React from "react"
import { Button, CardBody, CardFooter, Form, FormGroup, Input, Label} from 'reactstrap';
import { Field, Formik } from 'formik';
import * as Yup from 'yup';
import { FieldErrorMessage } from '../../../genui';

class ParameterField extends React.Component {
  CTYPE_TO_COMPONENT = {
    string: name => <Field name={name} as={Input} type="text"/>,
    integer: name => <Field name={name} as={Input} type="number"/>,
    float: name => <Field name={name} as={Input} type="number" step="0.01"/>
  };

  render() {
    return this.CTYPE_TO_COMPONENT[this.props.parameter.contentType](this.props.parameter.name)
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
                    <Label htmlFor="name">Name</Label>
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
                    </FormGroup>
                    : null
                  }

                  {
                    this.props.parameters.map(param => (
                      <FormGroup key={param.name}>
                        <Label htmlFor={param.name}>{param.name}</Label>
                        <ParameterField parameter={param}/>
                        <FieldErrorMessage name={param.name}/>
                      </FormGroup>
                    ))
                  }
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
      name: 'New Model',
      description: 'Write more about this model if needed...',
      mode: this.modes[0].name,
    };
    const parameterDefaults = {};
    for (const param of this.parameters) {
      parameterDefaults[param.name] = this.CTYPE_TO_DEFAULT[param.contentType]
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
    };
    const parameterValidators = {};
    for (const param of this.parameters) {
      parameterValidators[param.name] = this.CTYPE_TO_VALIDATOR[param.contentType]
    }
    validationObj = Object.assign(
      validationObj,
      parameterValidators
    );
    const validationSchema = Yup.object(validationObj);

    return (
      <ModelForm
        initialValues={initialValues}
        validationSchema={validationSchema}
        modes={this.modes}
        parameters={this.parameters}
        handleCreate={this.props.handleCreate}
      />
    );
  }
}

export default ModelCreateForm;