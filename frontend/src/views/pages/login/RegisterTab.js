import { Button, Col, Container, Form, FormGroup, Input, Label, Row } from 'reactstrap';
import React from 'react';

export default function RegisterTab(props) {
  return (
    <Container>
      <h2>Register</h2>
      <Form>

        <Row>
          <Col md={6}>
            <FormGroup>
              <Label for="username">Username</Label>
              <Input type="test" name="username" id="username" placeholder="Username" />
            </FormGroup>

            <FormGroup>
              <Label for="displayName">Display Name</Label>
              <Input type="text" name="displayName" id="displayName" placeholder="Display Name" />
            </FormGroup>
          </Col>

          <Col md={6}>
            <FormGroup>
              <Label for="email">E-mail</Label>
              <Input type="email" name="email" id="email" placeholder="E-mail" />
            </FormGroup>

            <FormGroup>
              <Label for="password">Password</Label>
              <Input type="password" name="password" id="password" placeholder="password" />
            </FormGroup>
          </Col>
        </Row>

        <Button>Register</Button>
      </Form>
    </Container>
  )
}