import {
  Button,
  Card,
  CardBody, CardFooter,
  Col,
  Container,
  Form,
  FormGroup,
  Input,
  Label,
  Row,
  UncontrolledAlert,
} from 'reactstrap';
import React from 'react';
import { Field, Formik } from 'formik';
import * as Yup from 'yup';
import { FieldErrorMessage } from '../../../genui';
import { Redirect, useHistory } from 'react-router-dom';
import LogInManager from './LoginManager';

function UserInfo(props) {
  return props.user && (props.user.first_name || props.last_name) ? (
    <React.Fragment>{props.user.username} ({`${props.user.first_name} ${props.user.last_name}`})</React.Fragment>
  ) : (props.user ? <React.Fragment>{props.user.username}</React.Fragment> : null)
}

function LoginForm(props) {
  const history = useHistory();

  if (props.loginSuccess) {
    return <Redirect to={props.appPath}/>
  }

  return props.user ? (
      <Row>
        <Col sm={12}>
          <Card>
            <CardBody className="display-flex">
              <p>You are already signed in as: <strong className="text-success"><UserInfo {...props}/></strong>
              </p>
            </CardBody>
            <CardFooter>
              <Button color="primary" onClick={e => {
                e.preventDefault();
                const path = props.appPath;
                history.push(path);
              }}>Go to App</Button> <Button color="danger" onClick={e => {
                e.preventDefault();
                props.logoutUser(props.user)
              }}>Logout</Button>
            </CardFooter>
          </Card>
        </Col>
      </Row>
    ) : (
    <Container>
      <h2>Login</h2>
      {
        props.apiErrors.map((err, index) => <UncontrolledAlert key={index} color="danger">{err}</UncontrolledAlert>)
      }
      <Formik
        initialValues={{
          username: "",
          password: ""
        }}
        validationSchema={Yup.object().shape({
          username: Yup.string().required("Username is required."),
          password: Yup.string().required("Password is required."),
        })}
        onSubmit={props.onSubmit}
      >
        {
          formik => (
            <Form onSubmit={formik.handleSubmit}>
              <FormGroup>
                <Label htmlFor="username">Username</Label>
                <Field name="username" as={Input} type="text" placeholder="Username" />
              </FormGroup>
              <FieldErrorMessage name="username"/>

              <FormGroup>
                <Label htmlFor="password">Password</Label>
                <Field as={Input} type="password" name="password" placeholder="Password" />
              </FormGroup>
              <FieldErrorMessage name="password"/>

              <Button type="submit" color="primary" disabled={formik.isSubmitting}>Login</Button>
            </Form>
          )
        }
      </Formik>
    </Container>
  )
}

export default function LoginTab(props) {
  return (
    <Container>
      <LogInManager fetchUser={true} {...props} component={LoginForm} />
    </Container>
  )
}