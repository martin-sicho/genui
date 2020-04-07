import { Col, Row } from 'reactstrap';
import SelectedList from './SelectedList';
import React from 'react';

class SelectedListPage extends React.Component {

  shouldComponentUpdate(nextProps, nextState, nextContext) {
    return nextProps.selectedMolsRevision !== this.props.selectedMolsRevision;
  }

  render() {
    return (
      <React.Fragment>
        <h1>Selected Compounds</h1>
        <hr/>
        <Row>
          <Col sm={12}>
            {
              this.props.selectedMols.length > 0 ? <SelectedList {...this.props}/> : <p>Select compounds in the map to see details.</p>
            }
          </Col>
        </Row>
      </React.Fragment>
    )
  }
}

export default SelectedListPage