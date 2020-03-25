import React from 'react';
import { ComponentWithObjects } from '../../../../genui';
import MapSelect from './MapSelect';

function Maps(props) {
  const defaultMapClass = "Map";
  const molsetColorList = [
    // according to https://stackoverflow.com/questions/40673490/how-to-get-plotly-js-default-colors-list
    // if exhausted, the first color gets used again
    '#1f77b4',  // muted blue
    '#ff7f0e',  // safety orange
    '#2ca02c',  // cooked asparagus green
    '#d62728',  // brick red
    '#9467bd',  // muted purple
    '#8c564b',  // chestnut brown
    '#e377c2',  // raspberry yogurt pink
    '#7f7f7f',  // middle gray
    '#bcbd22',  // curry yellow-green
    '#17becf'   // blue-teal
  ];
  return (
    <ComponentWithObjects
      objectListURL={props.apiUrls.mapsRoot}
      emptyClassName={defaultMapClass}
      currentProject={props.currentProject}
      render={
        (mapObjects) => {
          const maps = mapObjects[defaultMapClass];
          return <MapSelect {...props} maps={maps} molsetColorList={molsetColorList}/>
        }
      }
    />
  )
}

export default Maps;
