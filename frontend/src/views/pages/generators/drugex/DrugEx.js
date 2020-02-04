import React from "react";

class DrugExPage extends React.Component {

  componentDidMount() {
    this.props.setPageTitle("DrugEx");
  }

  render() {
    return <div>This is the drugex page.</div>
  }
}

export default DrugExPage;