import React from 'react';
import { CardBody, CardHeader } from 'reactstrap';
import { CreateNewForm } from './CreateNewForm';

export class CreateNewCard extends React.Component {

  render() {
    return (
      <React.Fragment>
        <CardHeader>Create New Project</CardHeader>
        <CardBody className="scrollable">
          <CreateNewForm {...this.props}/>
        </CardBody>
      </React.Fragment>
    );
  }
}