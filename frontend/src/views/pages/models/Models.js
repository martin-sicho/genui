import React from "react";
import { ComponentWithObjects, ComponentWithResources, ModelsPage } from '../../../genui';
import ModelCard from './ModelCard';
import ModelCardNew from './ModelCardNew';

function Models(props) {
  const resources = {
    algorithmChoices : new URL('algorithms/', props.apiUrls.qsarRoot),
    descriptors: new URL('metrics/', props.apiUrls.qsarRoot),
    metrics: new URL('descriptors/', props.apiUrls.qsarRoot)
  };
  return (
    <ComponentWithResources definition={resources}>
      {
        (allLoaded, resources) => (
          allLoaded ? <ComponentWithObjects
            objectListURL={new URL('all/', props.apiUrls.compoundSetsRoot)}
            {...props}
            render={
              (
                ...args
              ) => {
                const [compoundSets] = [...args];
                return (<ModelsPage
                  {...props}
                  {...resources}
                  modelClass="QSARModel"
                  listURL={new URL(`models/`, props.apiUrls.qsarRoot)}
                  modelComponent={ModelCard}
                  newModelComponent={ModelCardNew}
                  compoundSets={compoundSets}
                />)
              }
            }
          /> : <div>Loading...</div>
        )
      }
    </ComponentWithResources>
  );
}

export default Models;