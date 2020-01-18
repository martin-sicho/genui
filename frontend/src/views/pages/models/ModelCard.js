import React from "react";
import { CardBody, CardHeader } from 'reactstrap';

class ModelCard extends React.Component {

  render() {
    return (
      <React.Fragment>
        <CardHeader>Some Model</CardHeader>
        <CardBody>Data</CardBody>
      </React.Fragment>
    )
  }
}

export default ModelCard;