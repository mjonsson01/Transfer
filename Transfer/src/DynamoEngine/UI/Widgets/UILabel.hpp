// File: Transfer/src/DynamoEngine/UI/Widgets/UILabel.hpp

#pragma once

// Custom Imports
#include "DynamoEngine/UI/UIElement.hpp"

// Standard Library Imports
#include <functional>
#include <string>

namespace DynamoEngine
{
class UILabel : public UIElement
{
  public:
    using TextSource = std::function<std::string()>;

    explicit UILabel(std::string text) : m_text(std::move(text)) {}

    void setText(std::string text) { m_text = std::move(text); }
    void setTextSource(TextSource source) { m_text_source = std::move(source); }
    const std::string& text() const { return m_text; }

    void update(float delta_seconds) override;
    void draw(UIGeometryBuilder& builder) const override;

  private:
    std::string m_text;
    TextSource m_text_source;
};
} // namespace DynamoEngine
