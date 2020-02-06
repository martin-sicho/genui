import React from "react";
import { ModelsPage, ComponentWithObjects } from '../../../../genui';
import { DrugExAgentCreateCard, DrugExNetCreateCard } from './ModelCreateCards';
import { DrugExAgentCard, DrugExNetCard } from './ModelCards';

class DrugExModelList extends React.Component {

  render() {
    const modelClass = this.props.modelClass;
    const selectedToAdd = this.props.algorithmChoices.find(element => element.name === modelClass);

    return (
      <React.Fragment>
        <h1>{this.props.title}</h1>
        <hr/>
        <ModelsPage
          {...this.props}
          headerComponent={null}
          selectedToAdd={selectedToAdd}
          newCardSetup={this.props.newCardSetup}
          cardSetup={this.props.cardSetup}
        />
      </React.Fragment>
    )
  }
}

class DrugExPage extends React.Component {

  constructor(props) {
    super(props);

    this.netsUrl = new URL('networks/', this.props.apiUrls.drugexRoot);
    this.agentsUrl = new URL('agents/', this.props.apiUrls.drugexRoot);
    this.qsarUrl = new URL(`models/`, this.props.apiUrls.qsarRoot);
    this.INIT_MAP = {
      DrugExNetwork : {
        listURL : this.netsUrl,
        newModelComponent: DrugExNetCreateCard,
        modelComponent: DrugExNetCard,
        title: "DrugExNetworks"
      },
      DrugExAgent: {
        listURL : this.agentsUrl,
        newModelComponent: DrugExAgentCreateCard,
        modelComponent: DrugExAgentCard,
        title: "DrugExAgents"
      },
    }
  }

  componentDidMount() {
    this.props.setPageTitle("DrugEx");
  }

  render() {
    return (<div className="drugex-models-grids">
        {
          Object.keys(this.INIT_MAP).map(ModelClass => {
            const data = this.INIT_MAP[ModelClass];

            if (ModelClass === "DrugExAgent") {
              return (
                <div key={ModelClass} className={ModelClass}>
                  <ComponentWithObjects
                    objectListURL={this.qsarUrl}
                    emptyClassName={"QSARModel"}
                    {...this.props}
                    render={
                      (QSARModels) => {
                        const models = QSARModels["QSARModel"];
                        return (models.length > 0 ? <DrugExModelList
                            {...this.props}
                            {...data}
                            netsUrl={this.netsUrl}
                            agentsUrl={this.agentsUrl}
                            qsarUrl={this.qsarUrl}
                            modelClass={ModelClass}
                            environments={models}
                            newCardSetup={{
                              h : {"md" : 14, "sm" : 14},
                              w : {"md" : 1, "sm" : 1},
                              minH : {"md" : 3, "sm" : 3},
                            }}
                            cardSetup={{
                              h : {"md" : 12, "sm" : 12},
                              w : {"md" : 1, "sm" : 1},
                              minH : {"md" : 3, "sm" : 3},
                            }}
                          /> : <div>Loading...</div>
                        )
                      }
                    }
                  />
                </div>
              )
            }

            return (
              <div key={ModelClass} className={ModelClass}>
                <DrugExModelList
                  {...this.props}
                  {...data}
                  netsUrl={this.netsUrl}
                  agentsUrl={this.agentsUrl}
                  qsarUrl={this.qsarUrl}
                  modelClass={ModelClass}
                  newCardSetup={{
                    h : {"md" : 12, "sm" : 12},
                    w : {"md" : 1, "sm" : 1},
                    minH : {"md" : 3, "sm" : 3},
                  }}
                  cardSetup={{
                    h : {"md" : 11, "sm" : 11},
                    w : {"md" : 1, "sm" : 1},
                    minH : {"md" : 3, "sm" : 3},
                  }}
                />
              </div>
              )
          })
        }
      </div>
    )
  }
}

export default DrugExPage;