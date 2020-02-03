import React from "react";
import { ComponentWithObjects, ComponentWithResources, ModelCardNew, ModelsPage } from '../../../genui';
import ModelCard from './ModelCard';
import ModelCreateForm from './CreateForm';

function NewQSARCard (props) {
  let molsets = [];
  Object.keys(props.compoundSets).forEach(
    (key) => molsets = molsets.concat(props.compoundSets[key])
  );

  return <ModelCardNew {...props} molsets={molsets} formComponent={ModelCreateForm}/>
}

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
                const compoundSetsAvailable = !(Object.keys(compoundSets).length === 0 && compoundSets.constructor === Object);
                return (compoundSetsAvailable ? <ModelsPage
                  {...props}
                  {...resources}
                  modelClass="QSARModel"
                  listURL={new URL(`models/`, props.apiUrls.qsarRoot)}
                  modelComponent={ModelCard}
                  newModelComponent={NewQSARCard}
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