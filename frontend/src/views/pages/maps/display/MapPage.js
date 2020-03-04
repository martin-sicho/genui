import {
  Card,
  CardBody,
  Col,
  DropdownItem,
  DropdownMenu,
  DropdownToggle,
  Row,
  UncontrolledDropdown,
} from 'reactstrap';
import React from 'react';
import Map from './Map';
import {MolsByMolsets, MolsByMolsetsTabs } from './MolsByMolsets';
import { MoleculeDetail } from '../../../../genui';
import {ComponentWithResources} from '../../../../genui/';

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

function MolActivityDetail(props) {
  return (
    <React.Fragment>
      <h4>Activities</h4>
      Some summary of the activity data for this molecule.
    </React.Fragment>
  )
}

class MapSelectorPage extends React.Component {

  constructor(props) {
    super(props);

    this.state = {
      selectedMap : null,
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

  render() {
    const selectedMap = this.state.selectedMap ? this.state.selectedMap : this.props.maps[0];
    const molsets = selectedMap.molsets;
    const definition = {};
    molsets.forEach(molset => {
      molset.activities.forEach(activitySetID => {
        definition[activitySetID] = new URL(`${activitySetID}/`, this.props.apiUrls.activitySetsRoot)
      })
    });
    const Comp = this.props.component;
    return (
      <ComponentWithResources
        definition={definition}
      >
        {
          (allLoaded, activitySets) => {
            return allLoaded ? (
              <Comp
                {...this.props}
                activitySets={activitySets}
                molsets={molsets}
                map={selectedMap}
              />
            ) : <div>Loading...</div>
          }
        }
      </ComponentWithResources>
    );
  }
}

class MapHandlers extends React.Component {
  constructor(props) {
    super(props);

    this.state = {
      selectedMols : [],
      selectedPoints : [],
      hoverMol: null
    }
  }

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
    const Comp = this.props.component;
    return (<Comp
      {...this.props}
      {...this.state}
      onMolsSelect={this.handleMolsSelect}
      onMolHover={this.handleMolHover}
    />)
  }
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