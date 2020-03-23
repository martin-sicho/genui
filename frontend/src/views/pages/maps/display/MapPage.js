import {
  Card,
  CardBody,
  Col,
  Row,
} from 'reactstrap';
import React from 'react';
import Map from './Map';
import { CompoundOverview } from '../../../../genui';

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

  handleDeselect = () => {
    if (this.props.onMolsSelect) {
      this.props.onMolsSelect([]);
    }
    if (this.props.onPointsSelect) {
      this.props.onPointsSelect([]);
    }
  };

  handleMolHover = (mol, point) => {
    this.setState({
      hoverMol : mol
    })
  };

  render() {
    let hoverMol = this.state.hoverMol;
    const selectedMap = this.props.selectedMap;
    const selectedMols = this.props.selectedMols;

    if (selectedMols.length === 1) {
      hoverMol = selectedMols[0];
    }
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
                  onDeselect={this.handleDeselect}
                  onMolHover={this.handleMolHover}
                />
              </CardBody>
            </Card>
          </Col>

          <Col md={3} sm={2}>
            {
              hoverMol ? (
                <CompoundOverview
                  {...this.props}
                  map={selectedMap}
                  mol={hoverMol}
                />
              ) : <div><p>Hover over a point in the map to see more.</p></div>
            }
          </Col>
        </Row>
      ) : <div>Select a map to display from the menu.</div>
    )
  }
}

export default MapPage;