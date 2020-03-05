import React from 'react';
import { ComponentWithResources } from '../../../../genui';
import { DropdownItem, DropdownMenu, DropdownToggle, UncontrolledDropdown } from 'reactstrap';

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

export default MapSelectorPage;