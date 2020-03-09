import React from 'react';
import { ComponentWithObjects } from '../../../../genui';
import MapSelect from './MapSelect';

function Maps(props) {
  const defaultMapClass = "Map";
  return (
    <ComponentWithObjects
      objectListURL={props.apiUrls.mapsRoot}
      emptyClassName={defaultMapClass}
      currentProject={props.currentProject}
      render={
        (mapObjects) => {
          const maps = mapObjects[defaultMapClass];
          return <MapSelect {...props} maps={maps}/>
        }
      }
    />
  )
}

export default Maps;
