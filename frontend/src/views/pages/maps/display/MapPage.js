import {
  Card,
  CardBody,
  Col,
  Row,
} from 'reactstrap';
import React from 'react';
import Map from './Map';
import { MoleculeDetail } from '../../../../genui';
import MapHandlers from './MapHandlers';
import MapSelectorPage from './MapSelectorPage';

class MapsPageComponents extends React.Component {

  render() {
    const hoverMol = this.props.hoverMol;
    const selectedMap = this.props.selectedMap;
    return (
      selectedMap ? (
        <Row>
          <Col md={9} sm={10}>
            <Card>
              <CardBody>
                <Map
                  {...this.props}
                  map={selectedMap}
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
                        map={selectedMap}
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
      ) : <div>Select a map to display from the menu.</div>
    )
  }
}

function MapPage(props) {

  const PageImpl = props => {
    return <MapHandlers {...props} component={MapsPageComponents}/>
  };
  return (
    <MapSelectorPage {...props} component={PageImpl}/>
  )
}

export default MapPage;