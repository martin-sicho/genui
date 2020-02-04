import React from "react";
import {
  ComponentWithObjects,
  ComponentWithResources,
  ModelsPage,
} from '../../../genui';
import QSARModelCard from './ModelCard';
import QSARModelCreateCard from './QSARModelCreateCard';

function Models(props) {
  const resources = {
    algorithmChoices : new URL('algorithms/', props.apiUrls.qsarRoot),
    descriptors: new URL('descriptors/', props.apiUrls.qsarRoot),
    metrics: new URL('metrics/', props.apiUrls.qsarRoot)
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
                const compoundSetsAvailable = !(Object.keys(compoundSets).length === 0 && compoundSets.constructor === Object);
                return (compoundSetsAvailable ? <ModelsPage
                  {...props}
                  {...resources}
                  modelClass="QSARModel"
                  listURL={new URL(`models/`, props.apiUrls.qsarRoot)}
                  modelComponent={QSARModelCard}
                  newModelComponent={QSARModelCreateCard}
                  compoundSets={compoundSets}
                /> : <div><p>There are currently no compound sets. You need to create one before building a QSAR model.</p></div>)
              }
            }
          /> : <div>Loading...</div>
        )
      }
    </ComponentWithResources>
  );
}

export default Models;