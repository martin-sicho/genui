import {
  Card,
  CardBody,
  Col,
  Row,
} from 'reactstrap';
import React from 'react';
import Map from './Map';
import { MoleculeDetail } from '../../../../genui';

class MapPage extends React.Component {
  constructor(props) {
    super(props);

    this.state = {
      hoverMol: null
    }
  }

  handleMolsSelect = (mols, points) => {
    if (this.props.onMolsSelect) {
      this.props.onMolsSelect(mols);
    }
    if (this.props.onPointsSelect) {
      this.props.onPointsSelect(points);
    }
  };

  handleMolHover = (mol, point) => {
    this.setState({
      hoverMol : mol
    })
  };

  render() {
    const hoverMol = this.state.hoverMol;
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
                  onMolsSelect={this.handleMolsSelect}
                  onMolHover={this.handleMolHover}
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

export default MapPage;