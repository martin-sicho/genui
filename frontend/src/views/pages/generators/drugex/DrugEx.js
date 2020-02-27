import React from "react";
import { ComponentWithObjects, ComponentWithResources } from '../../../../genui';
import DrugExPage from './DrugExPage';

class DrugEx extends React.Component {

  render() {
    const resources = {
      algorithmChoices : new URL('algorithms/', this.props.apiUrls.generatorsRoot),
      metrics: new URL('metrics/', this.props.apiUrls.generatorsRoot)
    };

    return (
      <ComponentWithResources definition={resources}>
        {
          (allLoaded, resources) => (
            allLoaded ? <ComponentWithObjects
              objectListURL={new URL('all/', this.props.apiUrls.compoundSetsRoot)}
              emptyClassName={"MolSet"}
              {...this.props}
              render={
                (
                  compoundSets
                ) =>
                  (<DrugExPage
                    {...this.props}
                    {...resources}
                    compoundSets={compoundSets}
                  />)
              }
            /> : <div>Loading...</div>
          )
        }
      </ComponentWithResources>
    )
  }
}

export default DrugEx;