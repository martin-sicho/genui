import {
  Card,
  CardBody,
  Col,
  Row,
} from 'reactstrap';
import React from 'react';
import Map from './Map';
import { MolsByMolsets, MolsByMolsetsTabs } from './MolsByMolsets';
import { MoleculeDetail } from '../../../../genui';
import MapHandlers from './MapHandlers';
import MapSelectorPage from './MapSelectorPage';

function MolActivityDetail(props) {
  return (
    <React.Fragment>
      <h4>Activities</h4>
      Some summary of the activity data for this molecule.
    </React.Fragment>
  )
}

class MapsPageComponents extends React.Component {

  render() {
    const hoverMol = this.props.hoverMol;
    const selectedMap = this.props.map;
    return (
      selectedMap ? (

        <React.Fragment>
          <Row>
            <Col md={8} sm={10}>
              <Card>
                <CardBody>
                  <Map
                    {...this.props}
                  />
                </CardBody>
              </Card>

            </Col>

            <Col md={4} sm={2}>
              {
                hoverMol ? (
                  <Row>
                    <Col md={6} sm={4}>
                      <MoleculeDetail
                        {...this.props}
                        mol={hoverMol}
                      />
                    </Col>
                    <Col md={6} sm={8}>
                      <MolActivityDetail
                        {...this.props}
                      />
                    </Col>
                  </Row>
                ) : null
              }

              <hr/>

              <div>Activity and summary stats for the displayed molecule sets...</div>
            </Col>
          </Row>
          <Row>
            <Col sm={12}>
              <h1>Selected Compounds</h1>
              <hr/>
              <MolsByMolsets
                {...this.props}
                component={MolsByMolsetsTabs}
              />
            </Col>
          </Row>
        </React.Fragment>

      ) : <div>Select a map to display from the menu.</div>
    )
  }
}

function MapsPage(props) {

  const PageImpl = props => {
    return <MapHandlers {...props} component={MapsPageComponents}/>
  };
  return (
    <MapSelectorPage {...props} component={PageImpl}/>
    )
}

export default MapsPage;