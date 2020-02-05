import React from "react";
import { ComponentWithResources } from '../../../../genui';
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
            allLoaded ? <DrugExPage
              {...this.props}
              {...resources}
            /> : <div>Loading...</div>
          )
        }
      </ComponentWithResources>
    )

    // return (
    //   <ComponentWithObjects
    //     {...this.props}
    //     objectListURL={new URL('networks/', this.props.apiUrls.drugexRoot)}
    //     emptyClassName="DrugExNet"
    //     render={
    //       (drExNetworks, createDrExNetworkList, createDrExNetwork, deleteDrExNetwork, updateDrExNetwork) => {
    //         return (
    //           <ComponentWithObjects
    //             {...this.props}
    //             objectListURL={}
    //             emptyClassName="DrugExAgent"
    //             render={
    //               (drExAgents, createDrExAgentList, createDrExAgent, deleteDrExAgent, updateDrExAgent) => {
    //                 return (
    //                   <DrugExPage
    //                     {...this.props}
    //                     models={
    //                       Object.assign(drExNetworks, drExAgents)
    //                     }
    //                     createDrExNetList={createDrExNetworkList}
    //                     createDrExNet={createDrExNetwork}
    //                     deleteDrExNet={deleteDrExNetwork}
    //                     updateDrExNet={updateDrExNetwork}
    //                     createDrExAgentList={createDrExAgentList}
    //                     createDrExAgent={createDrExAgent}
    //                     deleteDrExAgent={deleteDrExAgent}
    //                     updateDrExAgent={updateDrExAgent}
    //                   />
    //                 )
    //               }
    //             }
    //           />
    //         )
    //       }
    //     }
    //   />
    // )
  }
}

export default DrugEx;