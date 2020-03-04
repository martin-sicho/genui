import React from 'react';
import { ComponentWithObjects } from '../../../../genui';
import MapsPage from './MapPage';

const MapDashboard = (props) => {
  const defaultMapClass= "Map";
  return (
    <ComponentWithObjects
      objectListURL={props.apiUrls.mapsRoot}
      emptyClassName={defaultMapClass}
      currentProject={props.currentProject}
      render={
        (mapObjects) => {
          const maps = mapObjects[defaultMapClass];
          return maps.length > 0 ? (
            <MapsPage {...props} maps={maps}/>
          ) : <div>Loading...</div>
        }
      }
    />
  );
};

export default MapDashboard;
