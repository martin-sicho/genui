import React from 'react';
import { ComponentWithObjects, ModelGrid } from '../../index';

class ModelsPage extends React.Component {

  constructor(props) {
    super(props);

    this.headerComponent = this.props.headerComponent;

    this.state = {
      selectedToAdd : null,
    }
  }

  handleAddNew = (model) => {
    this.setState({selectedToAdd : model})
  };

  componentDidMount() {
    const HeaderComp = this.headerComponent;
    if (HeaderComp) {
      this.props.onHeaderChange(<HeaderComp {...this.props} addChoices={this.props.algorithmChoices} onModelAdd={this.handleAddNew}/>);
    }
  }

  render() {
    const selectedToAdd = this.state.selectedToAdd ? this.state.selectedToAdd : this.props.selectedToAdd;

    return (
      <div className="models-page">
        <ComponentWithObjects
          {...this.props}
          emptyClassName={this.props.modelClass}
          objectListURL={this.props.listURL}
          render={
            (models, handleAddModelList, handleAddModel, handleModelDelete) => {
              return <ModelGrid
                {...this.props}
                models={models[this.props.modelClass]}
                chosenAlgorithm={selectedToAdd}
                handleAddModel={
                  (...args) => {
                    this.setState({selectedToAdd : this.props.selectedToAdd});
                    return handleAddModel(...args)
                  }
                }
                handleModelDelete={handleModelDelete}
              />
            }
          }
        />
      </div>
    );
  }
}

export default ModelsPage;