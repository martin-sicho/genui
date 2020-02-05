import React from "react";
import { ModelsPage } from '../../../../genui';
import { DrugExAgentCreateCard, DrugExNetCreateCard } from './ModelCreateCards';
import { DrugExAgentCard, DrugExNetCard } from './CreateCards';

class DrugExModelList extends React.Component {

  render() {
    return (
      <React.Fragment>
        <h1>{this.props.title}</h1>
        <hr/>
        <ModelsPage
          {...this.props}
          listURL={this.props.url}
          headerComponent={null}
          selectedToAdd={this.props.modelClass}
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
    this.INIT_MAP = {
      DrugExNet : {
        url : this.netsUrl,
        newModelComponent: DrugExNetCreateCard,
        modelComponent: DrugExNetCard,
        title: "DrugExNetworks"
      },
      DrugExAgent: {
        url : this.agentsUrl,
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
            return (
              <div key={ModelClass} className={ModelClass}>
                <DrugExModelList
                  {...this.props}
                  {...data}
                  modelClass={ModelClass}
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