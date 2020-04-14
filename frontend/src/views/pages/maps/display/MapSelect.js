import React from 'react';
import HeaderNav from './HeaderNav';
import { ComponentWithResources, IDsToResources } from '../../../../genui';
import MapTabs from './MapTabs';

class MapSelect extends React.Component {

  constructor(props) {
    super(props);

    this.state ={
      selectedMap: null
    }
  }

  componentDidUpdate(prevProps, prevState, snapshot) {
    const maps = this.props.maps;
    const selectedMap = this.state.selectedMap;
    if (maps.length > 0 && !selectedMap) {
      this.props.onHeaderChange(
        <HeaderNav
          {...this.props}
          maps={maps}
          onMapChoice={map => {
            this.setState({
              selectedMap: map,
            })
          }}/>
      );

      this.setState({
        selectedMap: maps[0]
      })
    }
  }

  render() {
    const maps = this.props.maps;
    const selected = this.state.selectedMap;

    if (selected) {
      const molsets = selected.molsets;
      if (molsets.length === 0) {
        return  <div>{selected.name} has no compounds associated with it. Maybe you deleted the compound sets associated with this map?</div>
      }

      const molsetsToColor = {};
      molsets.forEach((molset, index) => {
        if (index >= this.props.molsetColorList.length) {
          index = index % this.props.molsetColorList.length;
        }
        molsetsToColor[molset.id] = this.props.molsetColorList[index];
      });
      const resourcesDef = {};
      molsets.forEach(molset => {
        Object.assign(resourcesDef, IDsToResources(this.props.apiUrls.activitySetsRoot, molset.activities))
      });
      return (
        <ComponentWithResources
          selected={selected}
          definition={resourcesDef}
          updateCondition={
            (prevProps, nextProps) => prevProps.selected !== nextProps.selected
          }
        >
          {
            (allLoaded, activitySets) => {
              console.log(activitySets);
              return allLoaded ? (
                <MapTabs
                  {...this.props}
                  selectedMap={selected}
                  maps={maps}
                  molsets={molsets}
                  molsetsToColor={molsetsToColor}
                  activitySets={activitySets}
                />
              ) : <div>Loading...</div>
            }
          }
        </ComponentWithResources>
      )
    } else {
      return <div>Select a map from the menu above.</div>
    }
  }
}

export default MapSelect;