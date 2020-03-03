import { Card, CardBody, Col, DropdownItem, DropdownMenu, DropdownToggle, Row, UncontrolledDropdown } from 'reactstrap';
import React from 'react';
import Map from './Map';
import {MolsByMolsets, MolsByMolsetsTabs } from './MolsByMolsets';
import { DataPair, MoleculeData, MoleculeDetail } from '../../../../genui';

function HeaderNav(props) {
  return (<UncontrolledDropdown nav inNavbar>
    <DropdownToggle nav caret>
      Actions
    </DropdownToggle>
    <DropdownMenu right>
      <UncontrolledDropdown>
        <DropdownToggle nav>Show Map...</DropdownToggle>
        <DropdownMenu>
          {
            props.maps.map(choice =>
              (<DropdownItem
                key={choice.id}
                onClick={() => {props.onMapChoice(
                  choice
                )}}
              >
                {choice.name}
              </DropdownItem>)
            )
          }
        </DropdownMenu>
      </UncontrolledDropdown>
    </DropdownMenu>
  </UncontrolledDropdown>)
}

class MapsPage extends React.Component {

  constructor(props) {
    super(props);

    this.state = {
      selectedMap : null,
      selectedMols : [],
      selectedPoints : [],
      hoverMol: null
    }
  }

  componentDidMount() {
    this.props.onHeaderChange(
      <HeaderNav
        {...this.props}
        onMapChoice={this.handleMapChoice}/>
     );
  }

  handleMapChoice = (map) => {
    this.setState((prevState) => {
      return {
        selectedMap: map
      }
    })
  };

  handleMolsSelect = (mols, points) => {
    this.setState({
      selectedMols : mols,
      selectedPoints: points,
    })
  };

  handleMolHover = (mol, point) => {
    this.setState({
      hoverMol : mol
    })
  };

  render() {
    const selectedMap = this.state.selectedMap ? this.state.selectedMap : this.props.maps[0];
    const hoverMol = this.state.hoverMol;
    return (
      selectedMap ? (

        <React.Fragment>
          <Row>
            <Col md={8} sm={10}>
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

            <Col md={4} sm={2}>
              {
                hoverMol ? (
                  <Row>
                    <Col md={6} sm={4}>
                      <MoleculeDetail mol={hoverMol}/>
                    </Col>
                    <Col md={6} sm={8}>
                      <DataPair title="SMILES" data={hoverMol.smiles}/>
                      <DataPair title="InChiKey" data={hoverMol.inchiKey}/>
                      {
                        Object.keys(hoverMol.extraArgs).map(key => (
                          <DataPair key={key} title={key} data={hoverMol.extraArgs[key].toString()} />
                        ))
                      }
                    </Col>
                  </Row>
                ) : null
              }
              <div>Activity stats...</div>
            </Col>
          </Row>
          <Row>
            <Col sm={12}>
              <h1>Selected Compounds</h1>
              <hr/>
              <MolsByMolsets
                {...this.props}
                component={MolsByMolsetsTabs}
                map={selectedMap}
                selectedMols={this.state.selectedMols}
                selectedPoints={this.state.selectedPoints}
              />
            </Col>
          </Row>
        </React.Fragment>

      ) : <div>Select a map to display from the menu.</div>
    )
  }
}

export default MapsPage;