import React from "react";
import { ComponentWithObjects } from '../../../../genui';
import DrugExPage from './DrugExPage';

class DrugEx extends React.Component {

  render() {
    return (
      <ComponentWithObjects
        {...this.props}
        objectListURL={new URL('networks/', this.props.apiUrls.drugexRoot)}
        emptyClassName="DrugExNet"
        render={
          (drExNetworks, createDrExNetworkList, createDrExNetwork, deleteDrExNetwork, updateDrExNetwork) => {
            return (
              <ComponentWithObjects
                {...this.props}
                objectListURL={new URL('agents/', this.props.apiUrls.drugexRoot)}
                emptyClassName="DrugExAgent"
                render={
                  (drExAgents, createDrExAgentList, createDrExAgent, deleteDrExAgent, updateDrExAgent) => {
                    return (
                      <DrugExPage
                        {...this.props}
                        models={
                          Object.assign(drExNetworks, drExAgents)
                        }
                        createDrExNetList={createDrExNetworkList}
                        createDrExNet={createDrExNetwork}
                        deleteDrExNet={deleteDrExNetwork}
                        updateDrExNet={updateDrExNetwork}
                        createDrExAgentList={createDrExAgentList}
                        createDrExAgent={createDrExAgent}
                        deleteDrExAgent={deleteDrExAgent}
                        updateDrExAgent={updateDrExAgent}
                      />
                    )
                  }
                }
              />
            )
          }
        }
      />
    )
  }
}

export default DrugEx;