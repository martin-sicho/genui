import { Button, Col, Container, Form, FormGroup, Input, Label, Row, UncontrolledAlert } from 'reactstrap';
import React from 'react';
import { Formik, Field } from 'formik';
import * as Yup from 'yup';
import { FieldErrorMessage } from '../../../genui';

function RegisterForm(props) {
  if (props.registerSuccess) {
    return (
      <UncontrolledAlert color="success">
        Registration was successful.
        You can login with the credentials you provided now.
      </UncontrolledAlert>
    )
  }

  return (
    <Formik
      initialValues={{
        username: "",
        email: "",
        password1: "",
        password2: "",
      }}
      validationSchema={Yup.object().shape({
        username: Yup.string().required("Username is required."),
        email: Yup.string().email("E-mail is not in the correct format.").required("E-mail is required."),
        password1: Yup.string().required("Password is required."),
        password2: Yup.string().required("Password check is required."),
      })}
      onSubmit={props.onSubmit}
    >
      {
        formik => (
          <Form onSubmit={formik.handleSubmit}>
            <Row>
              <Col md={6}>
                <FormGroup>
                  <Label for="username">Username</Label>
                  <Field as={Input} type="test" name="username" placeholder="Username" />
                </FormGroup>
                <FieldErrorMessage name="username"/>

                <FormGroup>
                  <Label for="email">E-mail</Label>
                  <Field as={Input} type="email" name="email" placeholder="E-mail" />
                </FormGroup>
                <FieldErrorMessage name="email"/>
              </Col>

              <Col md={6}>
                <FormGroup>
                  <Label for="password1">Password</Label>
                  <Field as={Input} type="password" name="password1" placeholder="password" />
                </FormGroup>
                <FieldErrorMessage name="password1"/>

                <FormGroup>
                  <Label for="password2">Password Again</Label>
                  <Field as={Input} type="password" name="password2" placeholder="password (check)" />
                </FormGroup>
                <FieldErrorMessage name="password2"/>
              </Col>
            </Row>

            <Button type="Submit" color="primary" disabled={props.registerSuccess || formik.isSubmitting}>Register</Button>
          </Form>
        )
      }
    </Formik>
  )

}

class RegisterManager extends React.Component {

  constructor(props) {
    super(props);

    this.state = {
      registerSuccess: false,
    }
  }

  handleSubmit = (values, {setSubmitting, setErrors}) => {
    this.sendRegisterRequest(values, setErrors, setSubmitting);
  };

  sendRegisterRequest = (data, setErrors, setSubmitting) => {
    fetch(new URL('registration/', this.props.apiUrls.accountsRoot), {
      "headers": {
        "Accept": "application/json",
        "Content-Type": "application/json",
      },
      "body": JSON.stringify(data),
      "method": "POST"
    }).then(response => response.json())
      .then(data => {
        if (data.hasOwnProperty("key")) {
          this.setState({registerSuccess: true});
        } else {
          setErrors(data);
          setSubmitting(false);
        }
      })
      .catch(e => console.log(e));
  };

  render() {
    const Component = this.props.formComponent;
    return <Component {...this.props} {...this.state} onSubmit={this.handleSubmit}/>
  }
}

export default function RegisterTab(props) {
  return (
    <Container>
      <h2>Register</h2>
      <RegisterManager {...props} formComponent={RegisterForm}/>
    </Container>
  )
}