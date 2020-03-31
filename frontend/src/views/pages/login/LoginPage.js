import { Col, Container, Row } from 'reactstrap';
import { TabWidget } from '../../../genui';
import React from 'react';
import RegisterTab from './RegisterTab';
import LoginTab from './LoginTab';

export default function LoginPage(props) {

  const tabs = [
    {
      title: "Login",
      renderedComponent: LoginTab
    },
    {
      title: "Register",
      renderedComponent: RegisterTab
    }
  ];

  return (
    <Container>
      <h1>Sign In</h1>

      <hr/>

      <Row>
        <Col sm={12}>
          <TabWidget {...props} tabs={tabs}/>
        </Col>
      </Row>
    </Container>
  )
}