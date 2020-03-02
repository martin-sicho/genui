import { Card, CardBody, Col, DropdownItem, DropdownMenu, DropdownToggle, Row, UncontrolledDropdown } from 'reactstrap';
import React from 'react';
import Map from './Map';
import {MolsByMolsets, MolsByMolsetsTabs } from './MolsByMolsets';

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

  render() {
    const selectedMap = this.state.selectedMap ? this.state.selectedMap : this.props.maps[0];
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
                  />
                </CardBody>
              </Card>

            </Col>

            <Col md={4} sm={2}>
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