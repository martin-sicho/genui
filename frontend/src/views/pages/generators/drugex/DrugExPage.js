import React from "react";

class DrugExNetGrid extends React.Component {

  render() {
    return <div>This is the network grid.</div>
  }
}

class DrugExAgentGrid extends React.Component {

  render() {
    return <div>This is the agent grid</div>
  }
}

class DrugExPage extends React.Component {

  CLASS_TO_COMPONENT = {
    DrugExNet : DrugExNetGrid,
    DrugExAgent : DrugExAgentGrid
  };

  componentDidMount() {
    this.props.setPageTitle("DrugEx");
  }

  render() {
    const models = this.props.models;
    console.log(models);

    return (<div className="drugex-models-grids">
        {
          Object.keys(models).map(ModelClass => {
            if (this.CLASS_TO_COMPONENT.hasOwnProperty(ModelClass)) {
              const ModelComponent = this.CLASS_TO_COMPONENT[ModelClass];
              return (
                <div key={ModelClass} className={ModelClass}>
                  <ModelComponent {...this.props} models={models[ModelClass]} currentModelClass={ModelClass}/>
                </div>
              )
            } else {
              console.log(`Unknown class: ${ModelClass}`);
              return null;
            }
          })
        }
      </div>
    )
  }
}

export default DrugExPage;