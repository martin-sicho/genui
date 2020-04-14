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
      hoverMol: null,
      hoverOverview: null
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
    mol = !(this.props.selectedMols.length === 1) ? mol : this.props.selectedMols[0];
    if (!this.state.hoverMol || (mol.id !== this.state.hoverMol.id)) {
      this.setState({
        hoverMol : mol,
        hoverOverview: (props) => <CompoundOverview {...props} mol={mol} />
      })
    }
  };

  render() {
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
                  onDeselect={this.handleDeselect}
                  onMolHover={this.handleMolHover}
                />
              </CardBody>
            </Card>
          </Col>

          <Col md={3} sm={2}>
            {
              this.state.hoverOverview ? (
                <this.state.hoverOverview
                  {...this.props}
                  map={selectedMap}
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