import React from 'react';
import { Col, Row } from 'reactstrap';

class ChEMBLInfo extends React.Component {

  constructor(props) {
    super(props);

    this.molset = this.props.molset;
  }

  render() {
    return (
      <Row>
        <Col sm="12">
          <h4>Description</h4>
          <p>{this.molset.description}</p>
          <h4>Stats</h4>
          <p>Some stats...</p>
        </Col>
      </Row>
    );
  }
}

export default ChEMBLInfo;