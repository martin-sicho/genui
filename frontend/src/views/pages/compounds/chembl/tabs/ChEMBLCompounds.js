import React from 'react';
import { Col, Row } from 'reactstrap';
import { CompoundListFromAPI } from '../../../../../genui'

class ChEMBLCompounds extends React.Component {

  constructor(props) {
    super(props);

    this.molset = this.props.molset;
  }

  render() {
    return (
      <Row>
        <Col sm="12">
          <h4>Molecules in {this.molset.name}</h4>
          <CompoundListFromAPI {...this.props}/>
        </Col>
      </Row>
    );
  }
}

export default ChEMBLCompounds;