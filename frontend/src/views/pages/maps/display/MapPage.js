import { Col, DropdownItem, DropdownMenu, DropdownToggle, Row, UncontrolledDropdown } from 'reactstrap';
import React from 'react';
import Map from './Map';
import MapSidebar from './MapSidebar';

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
      selectedMap : null
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
    return (
      <Row>
        <Col md={8} sm={10}>
          <Map
            {...this.props}
            map={this.state.selectedMap}
          />
        </Col>

        <Col md={4} sm={2}>
          <MapSidebar
            {...this.props}
            currentMap={this.state.selectedMap}
          />
        </Col>
      </Row>
    )
  }
}

export default MapsPage;