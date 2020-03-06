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

class MapsPageComponents extends React.Component {

  render() {
    const hoverMol = this.props.hoverMol;
    const selectedMap = this.props.map;
    return (
      selectedMap ? (

        <React.Fragment>
          <Row>
            <Col md={9} sm={10}>
              <Card>
                <CardBody>
                  <Map
                    {...this.props}
                  />
                </CardBody>
              </Card>

            </Col>

            <Col md={3} sm={2}>
              {
                hoverMol ? (
                  <React.Fragment>
                    <Row>
                      <Col sm={12}>
                        <MoleculeDetail
                          {...this.props}
                          mol={hoverMol}
                        />
                      </Col>
                    </Row>
                    <hr/>
                    <Row>
                      <Col sm={12}>
                        <div>PhysChem props of the compound -> if more compounds selected remove the structure above and show some graphs!</div>
                      </Col>
                    </Row>
                  </React.Fragment>
                ) : null
              }
            </Col>
          </Row>
          <hr/>
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